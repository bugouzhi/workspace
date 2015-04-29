package org.Spectrums;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import rt.RTAligner;

import IO.GenericSpectrumReader;
import Utils.FileIOUtils;
import Utils.SpectrumFilter;

/**
 * Construct a spectral library from a set of PSMs
 * with RT alignment to a reference run
 * @author Jian
 *
 */
public class SpectrumLibConstructor {
	private String spectrumDir;
	private String PSMFile;
	private LargeHashMap resultMap;    //all PSMs results 
	private LargeHashMap rtAlignMap;   //store the RTalign model for each run
	private HashMap<Spectrum, String> fileMap;  //keep track of spectrum file the representative spectrum is from
	ResultColumnIndex index;
	private SpectrumLib referenceRun;
	private Object refRunKey;
	private ResultColumnIndex resultIndex;
	
	public SpectrumLibConstructor(String spectrumDir, String PsmFile, ResultColumnIndex indices){
		this.spectrumDir = spectrumDir;
		this.PSMFile = PsmFile;
		this.index = indices;
		this.resultIndex= indices;
		this.fileMap=new HashMap();
	}
	
	public SpectrumLib constructSpectralLibrary(String outfile, String logfile, SpectrumFilter filter, double tolerance){
		splitPSMsByFiles();
		this.referenceRun = getReferenceRun();
		buildLib();
		this.writeLib(outfile);
		this.getLibraryStats(logfile, filter, tolerance);
		return this.referenceRun;
	}
	
	
	/**
	 * Align RT between the reference and each run
	 * adding new peptide to the reference runs so there are points
	 * for alignment in subsequent runs
	 */
	private void buildLib(){
		NavigableMap<Integer, Object> orderedRuns = getIncludeOrder(this.referenceRun);
		this.rtAlignMap = new LargeHashMap(Utils.FileIOUtils.stripExtension(this.PSMFile+"_tmp_RT.map"));
		for(Iterator<Integer> it = orderedRuns.descendingKeySet().iterator(); it.hasNext();){
			Object curr = orderedRuns.get(it.next());
			SpectrumLib currLib = (SpectrumLib)this.resultMap.get(curr);
			Set intersect = Utils.SetUtils.getIntersect(this.referenceRun.getSpectrumLibrary().keySet(), 
					currLib.getSpectrumLibrary().keySet());
			System.out.println("reference size: " + this.referenceRun.getSpectrumLibrary().size());
			System.out.println("Intersection with reference: " + intersect.size());
			RTAligner align = new RTAligner(currLib, this.referenceRun);
			align.computeRegression();
//			System.out.println("Global mode");
//			align.getRTDiffInterval();
//			System.out.println("Local mode");
//			align.setPredictionMode(RTAligner.LOCAL);
//			align.getRTDiffInterval();
			this.rtAlignMap.put(curr, align);
			
			for(Iterator it2 = currLib.getSpectrumLibrary().keySet().iterator(); it2.hasNext();){
				String pep = (String)it2.next();
				List<Spectrum> specList = currLib.getSpectra(pep);
				Spectrum represent = this.getRepresentativeSpectrum(specList);
				for(int i = 0; i < specList.size(); i++){
					Spectrum s2 = specList.get(i);
					s2.rt = align.getAlignedRT(s2.rt);	
				}
				//adding new peptide and replaceing better representative spect in reference library 
				if(!intersect.contains(pep)){
					this.referenceRun.addSpectrum(represent);
					this.fileMap.put(represent, (String)curr);
				}else{
					List<Spectrum> currList = this.referenceRun.getSpectra(pep);
					Spectrum currRef = this.referenceRun.getSpectra(pep).get(0);
					if(currRef.score > represent.score){
						represent.rt = currRef.rt;    //we will use the reference RT instead of the predicted RT
						currList.remove(currRef);
						currList.add(represent);
						this.fileMap.remove(currRef);
						this.fileMap.put(represent, (String)curr);
					}
				}
				
			}
		}
	}
		
	
	/**
	 * Print out the library in MGF format
	 * @param outFile
	 */
	public  void writeLib(String outFile){
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		out.println(this.referenceRun);
		out.println();
		out.flush();
		out.close();
	}
	
	/**
	 * Get library stat that evaluate quality of libraries
	 * @param statFile
	 * @param filter
	 * @param tolerance
	 */
	public void getLibraryStats(String statFile, SpectrumFilter filter, double tolerance){
		List<Spectrum> spectList = this.referenceRun.getAllSpectrums();
		PrintStream statOut = Utils.FileIOUtils.getOutStream(statFile);
		statOut.print("#SpectrumFile\tScan#\tPrecursor\tCharge\tPeptide\tProtein\tRT\tScore\tMedianNoise\tSNR\tAnnotSig\tfractSig_Annotated\tfractAnnotated_isSig\tAnnotNoise\tunAnnotNoise\n");
		for(int i = 0; i < spectList.size(); i++){
			Spectrum s = spectList.get(i);
			String peptide = s.peptide;
			int charge = s.charge;
			Peptide p = new Peptide(peptide, charge);
			if(!this.fileMap.containsKey(s)){
				System.out.println("cannot find file");
			}
			String file = this.fileMap.get(s);
			statOut.print(file+"\t"+s.scanNumber + "\t" + s.parentMass + "\t" + s.charge + "\t" + s.peptide + "\t" + s.protein + "\t" + s.rt + "\t" + s.score+"\t");
			double[] stats = filter.computeAnnotatedSNR(s, tolerance);
			for(int j = 0; j < stats.length; j++){
				statOut.print(stats[j] +"\t");
			}
			statOut.println();
			if(Math.abs(p.getParentmass() - s.parentMass) > 1.2){
				//System.out.println(s);
				System.out.println("warning: parent mass not matching");
			}
		}
		statOut.flush();
		statOut.close();
	}
	
	
	
	/**
	 * Choose an order to gradually align other runs to the reference in RT space
	 * choosing run that share most peptides with reference so we can get
	 * most accurate alignment while gradually expanding our list of peptides in library
	 * @return
	 */
	private NavigableMap<Integer, Object> getIncludeOrder(SpectrumLib reference){
		NavigableMap<Integer, Object> orderedRuns = new TreeMap<Integer, Object>();
		for(Iterator<Object> it = this.resultMap.getKeys().iterator(); it.hasNext();){
			Object key = it.next();
			SpectrumLib lib = (SpectrumLib)this.resultMap.get(key);
			//System.out.println("currentlib size: " + lib.getSpectrumLibrary().keySet().size());
			//if(lib != reference){
			if(key != this.refRunKey){
				Set intersect =Utils.SetUtils.getIntersect(reference.getSpectrumLibrary().keySet(), 
						lib.getSpectrumLibrary().keySet());
				System.out.println("result " + key + " overlap with reference: " 
						+ intersect.size());
				orderedRuns.put(intersect.size(), key);
			}
		}
		return orderedRuns;
	}
	
	
	/**
	 * choose the reference run for RT-alignment to construct library
	 * @return the run which is defined as the run with most identified peptides
	 * the reference run contain a representative spectrum for each peptide
	 */
	private SpectrumLib getReferenceRun(){
		int maxNumPeps = 0;
		SpectrumLib reference = null;
		Object refKey = null;
		for(Iterator<Object> it = this.resultMap.getKeys().iterator(); it.hasNext();){
			Object key = it.next();
			SpectrumLib lib = (SpectrumLib)this.resultMap.get(key);
			System.out.println("result from " + key + "\thas peptides: " + lib.getSpectrumLibrary().keySet().size());
			int numPeps = lib.getSpectrumLibrary().keySet().size();
			reference = numPeps > maxNumPeps ? lib : reference;
			refKey = numPeps > maxNumPeps ? key : refKey;
			maxNumPeps = numPeps > maxNumPeps ? numPeps : maxNumPeps;
		}
		SpectrumLib ret = new SpectrumLib();
		for(Iterator<String> it = reference.getSpectrumLibrary().keySet().iterator(); it.hasNext();){
			String pepKey = it.next();
			List<Spectrum> specs = reference.getSpectra(pepKey);
			Spectrum s = this.getRepresentativeSpectrum(specs);
			if(s!=null){
				//System.out.println("adding "  + s.peptide + "\t"  + s.score);
				ret.addSpectrum(s);
				this.fileMap.put(s, (String) refKey);
			}
		}
		this.refRunKey = refKey;
		//System.out.println("Using reference " + refKey + "with IDs " + maxNumPeps);
		return ret;
		
	}
	

	
	/**
	 * Given a list of spectra from the same peptide ion, choose one
	 * spectrum to be a representative
	 * @param specList
	 * @return the spectrum with lowest score since we assume score is some kind of probablitliy or pvalue from search tools
	 */
	private Spectrum getRepresentativeSpectrum(List<Spectrum> specList){
		Spectrum represent=null;
		double bestScore = this.resultIndex.getSortOrder() > 0  ? Double.MAX_VALUE : Double.MIN_VALUE;
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(this.resultIndex.getSortOrder() > 0){
				represent = s.score < bestScore ? s : represent;
				bestScore = s.score < bestScore ? s.score : bestScore;
			}else{
				represent = s.score > bestScore ? s : represent;
				bestScore = s.score > bestScore ? s.score : bestScore;
			}
		}
		return represent;
	}

	
	/**
	 * split all the PSMs by input file name (often corresponding to experimental run)
	 * store them in a table
	 */
	private void splitPSMsByFiles(){
		resultMap = new LargeHashMap(Utils.FileIOUtils.stripExtension(this.PSMFile+"_tmp.map"));
		List<Spectrum> specList = new ArrayList<Spectrum>();
		Map<String, GenericSpectrumReader> readers = new HashMap<String, GenericSpectrumReader>();
		try{
			BufferedReader buff = new BufferedReader(new FileReader(this.PSMFile));
			String line;
			String prevFile= "";
			//MZXMLReader reader = null;
			GenericSpectrumReader reader = null;
			int counter = 1; //one base index
			int pepInd = this.index.getPepInd();
			int chargeInd = this.index.getChargeInd();
			while((line=buff.readLine()) != null){
				//System.out.println("Processing: " + line);
				if(line.startsWith("#")){
					continue;
				}
				String[] tokens = line.split("\\t");
				String file = tokens[0];
				file = file.replaceAll("\\s+", "");
				file = FileIOUtils.getFileName(file);
				if(!file.equals(prevFile)){
					if(specList.size() > 0){
						SpectrumLib lib = new SpectrumLib(specList, true);
						resultMap.put(prevFile, lib);
						specList=new ArrayList<Spectrum>();
					}
					reader = new GenericSpectrumReader(spectrumDir+File.separator + file);
					prevFile = file;
				}
				Spectrum s = reader.readSpectrumByIndex(Integer.parseInt(tokens[this.index.getSpecInd()]));
				String peptide = tokens[pepInd];
				int charge = Integer.parseInt(tokens[chargeInd]);
				peptide = Utils.StringUtils.stripNeighborRes(peptide);
				Peptide p = new Peptide(peptide, charge);
				s.peptide=peptide;
				s.charge = charge;
				s.score = Double.parseDouble(tokens[this.index.getScoreInd()]);
				if(Math.abs(p.getParentmass() - s.parentMass) > 1.2){
					//System.out.println(s);
					System.out.println("warning: parent mass not matching");
				}
				s.peptide=peptide;
				s.charge = Integer.parseInt(tokens[6]);
				//s.scanNumber = counter++;
				specList.add(s);
				if(this.index.getProtInd() > 0) s.protein = tokens[this.index.getProtInd()];
				counter++;
				if(counter % 10000 == 0){
					System.out.println("Processed PSMs: " + counter);
				}
			}
			//putting last runs into table
			SpectrumLib lib = new SpectrumLib(specList, true);
			resultMap.put(prevFile, lib);
			System.out.println("Number of total runs: " + this.resultMap.getKeys().size());
			this.resultMap.finalize();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public static void testConstructLib(){
		String specDir = "f:/workspace/msdata/UPS_plus_background/UPS_Ecoli/IDA_combine/";
		String psmFile = "../mixture_linked/LibraryCreation/ACG_swathdevelopment_IDA_combined_1PSMFDR_msgfdb_sorted.txt";
		String libOutFile = "../mixture_linked/LibraryCreation/testLib.mgf";
		String logOutFile = "../mixture_linked/LibraryCreation/log.txt";
		SpectrumLibConstructor constructor = new SpectrumLibConstructor(specDir, psmFile, ResultColumnIndex.MSGFDB_INDEX);
		SpectrumFilter filter = new SpectrumFilter();
		SpectrumLib lib = constructor.constructSpectralLibrary(libOutFile, logOutFile, filter, 0.05);
	}
		
	public static void main(String[] args){
		testConstructLib();
		
	}
	
	
}
