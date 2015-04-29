package dia;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Vector;

import IO.MZXMLReader;
import UI.CommandLineParser;

import org.Spectrums.DecoySpectrumGenerator;
import org.Spectrums.Mass;
import org.Spectrums.ModifiedSpectrum;
import org.Spectrums.PTM;
import org.Spectrums.Peak;
import org.Spectrums.PeakIntensityComparator;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.SimpleProbabilisticScorer;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumLib;
import org.Spectrums.SpectrumMap;
import org.Spectrums.TheoreticalSpectrum;
import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;



/**
 * Initial try on searching SWATH data with M-SPLIT
 * @author Jian Wang	
 *
 */
public class SWATHMSPLITSearch {
	public static int CENTROID = 1;
	public static int PROFILE = 0;
	public double parent = 25; //Da
	public double fragment = 50; //ppm
	public int topPeakKept = 15;
	public double windowWidth = 25;
	public int topLibPeaks = 30;
	public double minCos = 0.705;
	public static double maxSelfCos = 0.708;
	public int SWATHCycle = 35;
	public int minMatchedPeaks = 8;
	public int neighborWin = 5;
	public int maxRTDiff = 42000000;
	int dataMode = SWATHMSPLITSearch.CENTROID;
	List<PTM[]> ptmList;
	List<double[]> swathwindows=new ArrayList<double[]>();
	public String queryFile;
	public String libraryFile;
	public String outFile;
	
	/**
	 * This is the main entry point to SWATHMSPLIT searches
	 * @param minScan
	 * @param maxScan
	 */
	public void startMSPLITSearch(int minScan, int maxScan){
		Iterator<Spectrum> reader = new MZXMLReader(queryFile, minScan);
		ConsensusSpectrumReader reader2 = new ConsensusSpectrumReader(queryFile);
//		ConsensusSpectrumReader reader2 = new ConsensusSpectrumReader((MZXMLReader)reader);
		reader2.cycle = this.SWATHCycle;
		reader2.tolerance = this.fragment / 1000; 
		BufferedWriter bo;
		//parameters for searching	
		Double libTolDa = this.fragment / 1000; //tolerance in Da for library spectrum
		try{		
		bo = new BufferedWriter(new FileWriter(outFile));
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		//generateSimpliedSpectLib(lib);
		lib.mergeSpectrum(libTolDa);
		lib.filterPeaks(topLibPeaks);
		lib.normIntensity();
		getLibRT(lib);
		//swapRT(lib);
		Iterator<Spectrum> iter;
		iter=reader;
		int counter = 0;
		long start = 0;
		boolean go = true;
		bo.write("#File\tScan#\tMz\tz\tPeptide\tMz\tz\tcosine\tName\t#Peak(Query)\t#Peaks(match)\t#shared\tfraction-matched\trelative-alpha\tIonCount");
		for(int i = 0; i < 20; i++){
			bo.write("\tstat"+i);
		}
		bo.write("\n");
		Spectrum MS1 = null;
		SpectrumMap MSMap = null;
		List<Spectrum> libSpectList = lib.getAllSpectrums();
		if(this.ptmList != null){
			//lib.removeModSpectra();
			libSpectList = lib.getAllSpectrums();
			ModifiedSpectrum.addModToLibrary(libSpectList, this.ptmList);
			//ModifiedSpectrum.addModToLibrary(libSpectList);
		}
		while(iter.hasNext()){
			Spectrum s = iter.next();
			if(s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}else{
				if(go){
					start = (new GregorianCalendar()).getTimeInMillis();
					go = false;
				}				
			}
			//System.out.println(s.scanNumber);
			if(dataMode == PROFILE){
			//	s.mergePeaks(s, 0.032);
			}
			
			s.windowFilterPeaks2(topPeakKept, windowWidth);
			double tolerance = libTolDa;
			s.mergePeaks(s, tolerance);
			s.sqrtSpectrum();
			//System.out.println("Spectrum is: " + s + "\n");
			TreeMap<Double,Spectrum> bestCands = bestPsimSpec(libSpectList, s, MSMap, minCos, parent, this.fragment, maxRTDiff, Mass.DIFF_PPM);
			if(this.swathwindows.size() > 0){
				double[] window = this.swathwindows.get(counter % this.swathwindows.size());
				//System.out.println("Precursor " + s.parentMass + "\twindow:\t" + window[0] + "\t" + window[1]);
				bestCands = bestPsimSpec(libSpectList, s, MSMap, minCos, window[0], window[1], this.fragment, maxRTDiff, Mass.DIFF_PPM);
			}else{
				bestCands = bestPsimSpec(libSpectList, s, MSMap, minCos, parent, this.fragment, maxRTDiff, Mass.DIFF_PPM);
			}
			Iterator it = bestCands.descendingKeySet().iterator();
			double maxInt = 0.0; //matches with maximum abundance
			while(it.hasNext()){
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				double projectInt = cand.projectedPeakIntensity(s, tolerance);
				maxInt = projectInt > maxInt ? projectInt : maxInt;
			}
						
			it = bestCands.descendingKeySet().iterator();
			boolean DEBUG = false;
			while(it.hasNext()){	
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				if(cand.protein == null){
					cand.protein = cand.spectrumName;
				}
				
				double sharePeaks = cand.sharePeaks(s, tolerance, DEBUG);
				if(sharePeaks >= minMatchedPeaks){
					double projectInt = cand.projectedPeakIntensity(s, tolerance);
					double projectInt2 = s.projectedPeakIntensity(cand, tolerance);
					double libInt = cand.magnitude();
					libInt = libInt*libInt;
					//System.out.println("name is: " + s.spectrumName);
					cand.filterPeaksByMass(180, 2500);
					List<Spectrum> neighs = reader2.getNeighborScans(s.scanNumber, neighborWin, neighborWin);
//					neighs.add(new Spectrum());
//					for(int i = 0; i < neighs.size(); i++){
//						Spectrum neigh = neighs.get(i);
//						neigh.mergePeaks(s, 0.032);
//						neigh.windowFilterPeaks2(topPeakKept, windowWidth);
//						neigh.mergePeaks(s, 0.05);
//					}
					double[] simTwoD = reader2.getProjectCosine(cand, neighs, fragment);		
					bo.write(queryFile + "\t" +  s.scanNumber + "\t" + s.parentMass +"\t" + s.charge +"\t" 
							+ cand.peptide + "\t"  + cand.parentMass + "\t" + cand.charge + "\t" + psim 
							+ "\t" + cand.protein + "\t" + s.getPeaks().size() + "\t" + cand.getPeak().size()
							+"\t" + sharePeaks + "\t" + projectInt2 + "\t" + projectInt/maxInt + "\t" + s.upperBound +"\t" + libInt
							+"\t" + s.rt + "\t" + cand.rt + "\t");
					
					//double projCosTop = libSpect.projectedCosine(query, 0.05);
					for(int i = 0; i < simTwoD.length; i++){
						bo.write(simTwoD[i] + "\t");
					}
					bo.write(psim*simTwoD[5] +"\t" + simTwoD[1]*simTwoD[5]);
					bo.write("\n");
					//System.out.println();
				}
			}
			counter++;
			if(counter % 1000 == 0){
				System.out.println("Finish searched: " + counter);
			}
		}
		bo.flush();
		bo.close();
		//System.out.println("start second pass " + s.scanNumber);	
			System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		}catch(IOException IOE){
			System.out.println(IOE.getMessage());
			IOE.printStackTrace();
		}

	}
	
	/**
	 * More for older legacy spectral library format
	 * where RT information store in headers.  This method
	 * populated the proper field in the spectrum object with RT information.
	 * Newer annotated MGF has this as a RTINSECOND field
	 * @param lib
	 */
	public static void getLibRT(SpectrumLib lib){
		List<Spectrum> specList = lib.getSpectrumList();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s1 = specList.get(i);
			if(s1.spectrumName.split("\\s+").length > 5 && s1.spectrumName.contains("Retention Time: ")){  //check if retention time is in lib
				//System.out.println("SpectruMName: " + s1.spectrumName + "\t" + s1.spectrumName.split("\\s+")[5]);
				String RT = s1.spectrumName.split("Retention Time: ")[1].split("\\s+")[0];
				//String RT = s1.spectrumName.split("Retention Time: ")[5];
				RT=RT.replaceAll("[PTS]", "");
				s1.rt = Double.parseDouble(RT);
			}	

		}
	}
	
	/**
	 * Swap the retention time of decoy library spectrum,
	 * so they do not stay at the same RT region as the target peptides
	 * we swap with another spectrum fall within the same SWATH window to keep
	 * the overall distribution of RT the same as target and the distribution of
	 * the candidate peptide per swath window roughly the same
	 * @param lib
	 */
	protected static void swapRT(SpectrumLib lib){
		List<Spectrum> specList = lib.getSpectrumList();
		TreeMap<Double, Spectrum> sortedList = new TreeMap<Double, Spectrum>();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			//s.parentMass = s.parentMass+Math.random()*65;
			//s.parentMass = s.parentMass-1*Math.random()*65;
			if(s.spectrumName.contains("DECOY")){
				sortedList.put(s.parentMass+Math.random(), s);
			}
		}
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.spectrumName.contains("DECOY")){
				s.rt = s.rt + 500;
//				Map<Double, Spectrum> submap = sortedList.subMap(s.parentMass - 12, s.parentMass + 12);
//				List<Spectrum> values = new ArrayList<Spectrum>();
//				values.addAll(submap.values());
//				int ind = (int)Math.floor(Math.random()*values.size());
//				if(submap.size() > 0){
//					Spectrum swap = values.get(ind);
//					double rt = s.rt;
//					s.rt = swap.rt;
//					swap.rt = rt;
//					System.out.println("swapping rt: "  + s.spectrumName + "\t" + swap.spectrumName);
//				}
			}
		}
	}
	
	/**
	 * Generated a simplified version of the library spectrum by keeping
	 * only b/y ions of the library spectrum, test its effect on identification
	 * @return
	 */
	protected static void generateSimpliedSpectLib(SpectrumLib lib){
		DecoySpectrumGenerator gen = new DecoySpectrumGenerator();
		String[] prefixIons = new String[]{"b"};//, "b-H20", "b-NH3"};
		String[] suffixIons = new String[]{"y"};//, "y-H20", "y-NH3"};
		List<Spectrum> specList = lib.getAllSpectrums();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			Spectrum simplified = gen.generateTarget(s, prefixIons, suffixIons);
			System.out.println("Full spect size: " + s.getPeak().size());
			s.setPeaks(simplified.getPeak());
			System.out.println("Simplified spect size: " + s.getPeak().size());
		}
		TheoreticalSpectrum.prefixIons = new String[]{"b", "b-H20", "b-NH3"};//, "b(iso)"};
		TheoreticalSpectrum.suffixIons = new String[]{"y", "y-H20", "y-NH3"};//, "y(iso)"};
		
	}
	
	/**
	 * Automatically learned the number of of scans per cycle from data
	 * It scans the spectrum file and find number of scans that separate spectrum
	 * with the same precursor mass, it assumes this is the swath cycle 
	 * @return
	 */
	private int getSWATHCycle(){
		MZXMLReader reader = new MZXMLReader(this.queryFile);
		int counter = 0;
		Spectrum first = null;
		while(reader.hasNext()){
			if(first == null){
				first = reader.next();
			}else{
				Spectrum next = reader.next();
				if(Math.abs(next.parentMass - first.parentMass) < 0.05){
					counter = next.scanNumber - first.scanNumber;
					break;
				}
			}
		}
		System.out.println("DIA cycle from file: " + counter);
		return counter;
	}
	
	protected static void testMSPLITSearch(int minScan, int maxScan, 
			double parentMassWindow, double fragmentMassTol, String queryFile, String libraryFile, String outFile){
		//queryFile = "../mixture_linked/msdata/UPS_Ecoli_Wiff/Duplicate_runs_201308/REP2/18482_REP3_500ng_Ecoli_NewStock2_SWATH_2.mzXML";
     	//libraryFile = "../mixture_linked/ACG_swathdevelopment_UPSEcoli_REP234_IDA_plusDecoy2.mgf";
		//outFile = "../mixture_linked/swath_msplit_out.txt";
		Iterator<Spectrum> reader = new MZXMLReader(queryFile);
		ConsensusSpectrumReader reader2 = new ConsensusSpectrumReader(queryFile);
		BufferedWriter bo;
		try{		
			bo = new BufferedWriter(new FileWriter(outFile));
		//reader.minInt=40;
		//reader.numNeighbors=2;
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		//lib.windowFilterPeaks(6, 50);
		//lib.removePeakByMass(10, 400);
		//lib.removeModSpectra();
		lib.filterPeaks(30);
		lib.mergeSpectrum(0.05);
		//SpectrumLib lib = new SpectrumLib(libraryFile, "splib");
		//SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		//self-test: iter = lib.getSpectrumList().iterator();
		lib.normIntensity();
		Iterator<Spectrum> iter;
		//Iterator<Spectrum> iter = lib.getAllSpectrums().iterator();
		//iter = new LargeSpectrumLibIterator(queryFile);
		iter=reader;
		int counter = 0;
		long start = 0;
		boolean go = true;
		bo.write("#File\tScan#\tMz\tz\tPeptide\tMz\tz\tcosine\tName\t#Peak(Query)\t#Peaks(match)\t#shared\tfraction-matched\trelative-alpha\tIonCount\n");
		Spectrum MS1 = null;
		SpectrumMap MSMap = null;
		List<Spectrum> libSpectList = lib.getAllSpectrums();
		//System.out.println("Unmod lib size: " + libSpectList.size());
		//ModifiedSpectrum.addModToLibrary(libSpectList);
		//System.out.println("Modded lib size: " + libSpectList.size());
		while(iter.hasNext()){
			Spectrum s = iter.next();
			minScan = 300;
			maxScan = 301150;
			if(s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}else{
				if(go){
					start = (new GregorianCalendar()).getTimeInMillis();
					go = false;
				}				
			}
			
//			int ms1Scan = ((MZXMLReader)reader).getPrevScan(s.scanNumber,1);
//			if(MS1 == null || MS1.scanNumber != ms1Scan){
//				MS1 = ((MZXMLReader)reader).getSpectrum(ms1Scan);
//				MS1.toRelIntensity();
//				MS1.filterPeaksByIntensity(0.01);
//				//System.out.println(MS1);	
//				MSMap = new SpectrumMap(MS1);
//			}

			//s.filterPeaksByIntensity(100);
			//System.out.println(s);
			//s.toRelIntensity();
			//s.filterPeaksByIntensity(30);
			//s.computePeaksZScore(0.5);
			//System.out.println(s);
			//s.filterPeaksByIntensity(50);
			//s.filterPeaksByRankScore(220);
			//s.filterPeaks(500);
			
			//s.deIsoPeaks(s, 0.03);
			s.windowFilterPeaks2(15, 25);
			//s.windowFilterPeaks2(2.0, 50);
			//s.filterPeaksByRankScore(3);
			//s.filterPeaksByIntensity(100);
			//s.filterPeaksByIntensity(50);
			s.mergePeaks(s, 0.05);
			s.sqrtSpectrum();
			//System.out.println(s);
			//System.out.println("Spectrum: " + s.scanNumber + "\t" + s.parentMass +"\t" + s.charge + "\t#peaks: " + s.getPeaks().size());
			//TreeMap<Double,Spectrum> bestCands = bestPsimSpec(lib.getAllSpectrums(), s, 1, 2, 0.05);
			//s.shiftSpectrumPPM(100);
			double tolerance = 0.05;
			TreeMap<Double,Spectrum> bestCands = bestPsimSpec(libSpectList, s, MSMap, 0.72, 25, tolerance);
			Iterator it = bestCands.descendingKeySet().iterator();
			double maxInt = 0.0; //matches with maximum abundance
			while(it.hasNext()){
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				double projectInt = cand.projectedPeakIntensity(s, tolerance);
				maxInt = projectInt > maxInt ? projectInt : maxInt;
			}
						
			it = bestCands.descendingKeySet().iterator();
			boolean DEBUG = false;
			while(it.hasNext()){	
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				if(cand.protein == null){
					cand.protein = cand.spectrumName;
				}
				
				double sharePeaks = cand.sharePeaks(s, tolerance, DEBUG);
				if(sharePeaks > 8){
					double projectInt = cand.projectedPeakIntensity(s, tolerance);
					double projectInt2 = s.projectedPeakIntensity(cand, tolerance);
					double libInt = cand.magnitude();
					libInt = libInt*libInt;
					//System.out.println("name is: " + s.spectrumName);
					String RT = "0";
					if(cand.spectrumName.split("\\s+").length > 5 && cand.spectrumName.contains("Reteintion Time")){  //check if retention time is in lib
						RT = cand.spectrumName.split("\\s+")[5];
						RT = RT.substring(2, RT.length()-3);
					}	
					double libRT = Double.parseDouble(RT);
					cand.filterPeaksByMass(180, 2500);
					List<Spectrum> neighs = reader2.getNeighborScans(s.scanNumber, 5, 5);
					double[] simTwoD = reader2.getProjectCosine(cand, neighs, 50);		
					bo.write(queryFile + "\t" +  s.scanNumber + "\t" + s.parentMass +"\t" + s.charge +"\t" 
							+ cand.peptide + "\t"  + cand.parentMass + "\t" + cand.charge + "\t" + psim 
							+ "\t" + cand.protein + "\t" + s.getPeaks().size() + "\t" + cand.getPeak().size()
							+"\t" + sharePeaks + "\t" + projectInt2 + "\t" + projectInt/maxInt + "\t" + s.upperBound +"\t" + libInt
							+"\t" + s.rt + "\t" + libRT + "\t");
					
					//double projCosTop = libSpect.projectedCosine(query, 0.05);
					for(int i = 0; i < simTwoD.length; i++){
						bo.write(simTwoD[i] + "\t");
					}
					bo.write(psim*simTwoD[5] +"\t" + simTwoD[1]*simTwoD[5]);
					bo.write("\n");
					//System.out.println();
				}
			}
		//System.out.println("start second pass " + s.scanNumber);	
		//=======================================second pass==================================================  		
			//it = bestCands.descendingKeySet().iterator();
			it = new ArrayList().iterator();
			while(it.hasNext()){	
				Double sim = (Double)it.next();
				Spectrum current = bestCands.get(sim);				
				double sharePeaks = current.sharePeaks(s, 0.1);
				if(sim > 0.8 && sharePeaks > 15){
					Spectrum best = bestCands.lastEntry().getValue();
					Spectrum s2 = s.removeSharePeaks(current, 0.05);
					TreeMap<Double,Spectrum> bestCands2 = bestPsimSpec(lib.getAllSpectrums(), s2, 0.7, 25, 0.05);
					it = bestCands2.descendingKeySet().iterator();
					maxInt = 0.0; //matches with maximum abundance
					while(it.hasNext()){
						Double psim = (Double)it.next();
						Spectrum cand = bestCands2.get(psim);
						double projectInt = cand.projectedPeakIntensity(s, 0.05);
						maxInt = projectInt > maxInt ? projectInt : maxInt;
					}
					
					
					it = bestCands2.descendingKeySet().iterator();
					while(it.hasNext()){	
						Double psim = (Double)it.next();
						Spectrum cand = bestCands2.get(psim);
						if(cand.protein == null){
							cand.protein = cand.spectrumName;
						}
						sharePeaks = cand.sharePeaks(s2, 0.1);
						double projectInt = cand.projectedPeakIntensity(s2, 0.05);
						double projectInt2 = s2.projectedPeakIntensity(cand, 0.05);
						if(sharePeaks > 10){
//							System.out.println(queryFile + "\t" +  s2.scanNumber + "\t" + s2.parentMass +"\t" + s2.charge +"\t" 
//									+ cand.peptide + "\t"  + cand.parentMass + "\t" + cand.charge + "\t" + psim 
//									+ "\t" + cand.protein + "\t" + s2.getPeaks().size() + "\t" + cand.getPeak().size() 
//									+"\t" + sharePeaks + "\t" + projectInt + "\t" + projectInt/maxInt + "\t" + s.upperBound +"\t" );
//							System.out.println(s);
//							System.out.println();
//							System.out.println(cand);
						}
					}

				}
				
			}			
			
			//TreeMap<Double,Spectrum> bestCands = bestPsimSpec(reader, s, 10, 25, 0.03);
//			System.out.println("best cand size: " + bestCands.size());
			//System.out.println("best score: " + bestCands.lastKey());

			counter++;
			if(counter % 1000 == 0){
				System.out.println("Finish searched: " + counter);
			}
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		}catch(IOException IOE){
			System.out.println(IOE.getMessage());
			IOE.printStackTrace();
		}

	}
	
	/**
	 * Use parrelel DDA run to help make ID
	 * RT information is used to restrict search, the best matching SWATH is considered match
	 */
	public static void targetedIdentification(){
		String diaFile = "../mixture_linked/msdata/UPS_Ecoli_Wiff/Duplicate_runs_201308/REP2/18488_REP3_40fmol_UPS1_1ug_Ecoli_NewStock2_SWATH_1.mzXML";
		String ddaFile = "../mixture_linked/msdata/UPS_Ecoli_Wiff/IDA_combine/18487_REP3_40fmol_UPS1_1ug_Ecoli_NewStock2_IDA_1.mzXML";
		String ddaIDs = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/IDA/UPSEcoli_REP3_newStock2_msgfdb_2_pep.txt";
		MZXMLReader reader = new MZXMLReader(diaFile);
		MZXMLReader reader2 = new MZXMLReader(ddaFile);
		List<String> results = Utils.FileIOUtils.createListFromFile(ddaIDs);
		//ConsensusSpectrumReader reader = new ConsensusSpectrumReader(queryFile);
		//reader.minInt=40;
		//reader.numNeighbors=3;
		Spectrum MS1 = null;
		SpectrumMap MSMap = null;
		SortedMap<Double, Spectrum> timeMap = new TreeMap<Double, Spectrum>();
		MSXMLParser parser = reader2.getParser();
		for(int i = 0; i < results.size(); i++){
			String result = results.get(i);
			if(result.startsWith("#")){
				continue;
			}
			String[] tokens = result.split("\\s+");
			int scan = Integer.parseInt(tokens[1]);
			Scan s = parser.rap(scan);
			if(s.getHeader().getMsLevel() == 2){
				Spectrum spect = reader2.getSpectrum(scan);
				spect.windowFilterPeaks2(6, 50);
				spect.filterPeaks(50);
				spect.mergePeaks(spect, 0.05);
				spect.sqrtSpectrum();
				spect.peptide = tokens[7];
				spect.charge = Integer.parseInt(tokens[6]);
				//System.out.println("mappped RT: " + spect.spectrumName + "\t" + spect.peptide + "\t" + SWATHUtils.getRT(s));
				timeMap.put(SWATHUtils.getRT(s), spect);
			}
		}
		System.out.println("built time map");
		Iterator<Spectrum> iter = reader;
		MSXMLParser parser2 = reader.getParser();
		int minScan = 100;
		int maxScan = 151000;
		while(iter.hasNext()){
			Spectrum s = iter.next();
			if(s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}
			
			int ms1Scan = ((MZXMLReader)reader).getPrevScan(s.scanNumber,1);
			if(MS1 == null || MS1.scanNumber != ms1Scan){
				MS1 = ((MZXMLReader)reader).getSpectrum(ms1Scan);
				MS1.toRelIntensity();
				MS1.filterPeaksByIntensity(0.01);
				//System.out.println(MS1);
				MSMap = new SpectrumMap(MS1);
			}

			Scan scan = parser2.rap(s.scanNumber);
			double time = SWATHUtils.getRT(scan);
			//System.out.println("query:  " + s.scanNumber + "\tret-tme:\t" + time);
			SortedMap<Double, Spectrum> sub = timeMap.subMap(new Double((time - 300.0)), 
					new Double((time + 300.0)));
			//System.out.println("sub size: " + sub.size());
			List<Spectrum> lib = new ArrayList<Spectrum>();
			for(Iterator<Double> it = sub.keySet().iterator(); it.hasNext();){
				Spectrum libEntry = timeMap.get(it.next());
				//System.out.println("libEntry is: " + libEntry.peptide);
				if(MSMap.checkChargedPeaks(libEntry.parentMass, 0.03, libEntry.charge) &&
						s.parentMass - 5 < libEntry.parentMass){
					
					//System.out.println("Adding libEntry: " + libEntry.peptide);
					lib.add(libEntry);
				}
			}
			//System.out.println("checked precursor");
			s.filterPeaks(1000);
			//s.filterPeaksByIntensity(30);
			//s.computePeaksZScore(0.6);
			//s.filterPeaksByRankScore(200);
			//s.windowFilterPeaks2(15, 25);
			s.mergePeaks(s, 0.05);
			s.sqrtSpectrum();
			double tolerance = 0.035;
			TreeMap<Double,Spectrum> bestCands = bestPsimSpec(lib, s, 0.1, 25, tolerance);
			Iterator it = bestCands.descendingKeySet().iterator();
			double maxInt = 0.0; //matches with maximum abundance
			while(it.hasNext()){
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				double projectInt = cand.projectedPeakIntensity(s, tolerance);
				maxInt = projectInt > maxInt ? projectInt : maxInt;
			}
			
			
			it = bestCands.descendingKeySet().iterator();
			while(it.hasNext()){	
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				if(cand.protein == null){
					cand.protein = cand.spectrumName;
				}
				double sharePeaks = cand.sharePeaks(s, tolerance);
				double projectInt = cand.projectedPeakIntensity(s, tolerance);
				double projectInt2 = s.projectedPeakIntensity(cand, tolerance);
				double libInt = cand.magnitude();
				libInt = libInt*libInt;
				if(sharePeaks > 5){
					double[] scores = localSim(cand, s, 5, tolerance);
					System.out.print(diaFile + "\t" +  s.scanNumber + "\t" + s.parentMass +"\t" + s.charge +"\t" 
							+ cand.peptide + "\t"  + cand.parentMass + "\t" + cand.charge + "\t" + psim 
							+ "\t" + cand.protein + "\t" + s.getPeaks().size() + "\t" + cand.getPeak().size()
							+"\t" + sharePeaks + "\t" + projectInt + "\t" + projectInt2 + "\t" + projectInt/maxInt + "\t" + s.upperBound +"\t" + libInt);
					for(int i = 0; i < scores.length; i++){
						System.out.print(scores[i] + "\t");
					}
					//System.out.println(s);
					System.out.println();
					//System.out.println(cand);
				}
			}
		}
	}
	
	
	protected static double[] localSim(Spectrum libSpect, Spectrum query, int topN, double tolerance){
		Spectrum copy = new Spectrum(libSpect);
		double[] scores = new double[topN+1];
		List<Peak> sortedPeakList = new Vector<Peak>();
		sortedPeakList.addAll(copy.getPeak());
		Collections.sort(sortedPeakList, PeakIntensityComparator.comparator);
		scores[0] = copy.projectedCosine(query, tolerance);
		for(int i = 1; i <= topN && copy.getPeak().size() > 0; i++){
			Peak p = sortedPeakList.get(sortedPeakList.size()-i);
			copy.getPeak().remove(p);
			scores[i] = copy.projectedCosine(query, tolerance);
		}
		return scores;
	}
	
	/**
	 * Search in reverse direction, for each peptide search the best spectrum that matches to it
	 * @param minScan
	 * @param maxScan
	 */
	public static void testReverseMSPLITSearch(int minScan, int maxScan){
		String queryFile = "../mixture_linked/msdata/UPS_Ecoli/40fmol_UPS1_1ugEcoli_SWATH.mzXML";
		//String libraryFile = "../mixture_linked/spectral_library/human_qtof_plusDecoy.sptxt";
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test2.mgf";
		Iterator<Spectrum> reader = new MZXMLReader(queryFile);
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		double parentMassTol = 25.0;
		//SpectrumLib lib = new SpectrumLib(libraryFile, "splib");
		//lib.normIntensity();
		Iterator<Spectrum> iter = reader;
		int counter = 0;
		long start = new GregorianCalendar().getTimeInMillis();
		Map<Spectrum, Spectrum> bestMap = new HashMap();
		Map<Spectrum, Double> scoreMap = new HashMap();
		List<Spectrum> specList = lib.getAllSpectrums();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s1 = specList.get(i);
			//s1.windowFilterPeaks(6, 25);
			s1.filterPeaks(30);
			s1.mergePeaks(s1, 0.05);
			s1.sqrtSpectrum();
			bestMap.put(s1, null);
			scoreMap.put(s1, -1000.0);
		}
		while(iter.hasNext()){
			Spectrum s = iter.next();
			if(s.scanNumber < minScan){
				continue;
			}
			if(s.scanNumber > maxScan){
				break;
			}
			s.windowFilterPeaks2(15, 25);
			s.mergePeaks(s, 0.075);
			s.sqrtSpectrum();
			//System.out.println("Spectrum: " + s.scanNumber + "\t" + s.parentMass +"\t" + s.charge + "\t#peaks: " + s.getPeaks().size());
			double tolerance = 0.051;
			for(int i = 0; i < specList.size(); i++){
				Spectrum cand = specList.get(i);
				double massDelta = cand.parentMass - s.parentMass +4;
				if(massDelta > 0 && massDelta < 26){
					double sharePeaks = cand.sharePeaks(s, tolerance);
					double projectInt = cand.projectedPeakIntensity(s, tolerance);
					double projectInt2 = s.projectedPeakIntensity(cand, tolerance);
					double currScore = cand.projectedCosine(s, 0.05);
					//if(project)
					System.out.print(queryFile + "\t" +  s.scanNumber + "\t" + s.parentMass +"\t" + "\t" 
							+ cand.peptide + "\t"  + cand.parentMass + "\t" + cand.charge + "\t"  
							+ cand.protein + "\t" + s.getPeaks().size() + "\t" + cand.getPeak().size()
							+"\t" + sharePeaks + "\t" + projectInt + "\t" + projectInt2 +"\n");
					//if(scoreMap.get(s1) < currScore){
					//	scoreMap.put(s1, currScore);
					//	bestMap.put(s1, s);
					//}
					
				}
			}
			counter++;
			if(counter > 1000){
				//break;
			}
		}
		
		for(int i = 0; i < specList.size(); i++){
			Spectrum s1 = specList.get(i);
			Spectrum best = bestMap.get(s1);
			double score = scoreMap.get(s1);
			double sharePeaks = s1.sharePeaks(best, 0.1);
			double projectInt = s1.projectedPeakIntensity(best, 0.05);
			//System.out.println(s1.toString());
			if(best != null){
				//System.out.println(s1.spectrumName + "\t" + s1.peptide + "\t" + s1.parentMass + "\t" + s1.charge 
				//	+ "\t" +  best.scanNumber + "\t" + score + "\t" + best.getPeak().size() + "\t" + s1.getPeak().size() + "\t" + sharePeaks);
			}
		}

		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	/**
	 * Experimental method, this attempt to use theroretical spectrum to search for best swath spectrum match
	 * @param minScan
	 * @param maxScan
	 */
	protected static void testReverseTheoLibrarySearch(int minScan, int maxScan){
		String diaFile = "../mixture_linked/msdata/UPS_Ecoli//40fmol_UPS1_1ugEcoli_SWATH.mzXML";
		String ssmResult = "../mixture_linked/ACG_swathdevelopment_40fmol_UPS1_1ugEcoli_SWATH_vs_IDAlibrary_projcos_min0.6_top30_nodouble_pep.txt";
		MZXMLReader reader = new MZXMLReader(diaFile);
		List<String> results = Utils.FileIOUtils.createListFromFile(ssmResult);
		
		//ConsensusSpectrumReader reader = new ConsensusSpectrumReader(queryFile);
		//reader.minInt=40;
		//reader.numNeighbors=3;
		Spectrum MS1 = null;
		SpectrumMap MSMap = null;
		SortedMap<Double, Spectrum> timeMap = new TreeMap<Double, Spectrum>();
		MSXMLParser parser = reader.getParser();
		String singlescorer = "../mixture_linked/yeast_NIST_lib_singlepep_win10_25.o";
		RankBaseScoreLearner pcomp = RankBaseScoreLearner.loadComparator(singlescorer);
		SimpleProbabilisticScorer scorer3 = new SimpleProbabilisticScorer(pcomp);
		scorer3.setMatchTolerance(0.05);
		for(int i = 0; i < results.size(); i++){
			String result = results.get(i);
			String[] tokens = result.split("\\t");
			if(tokens.length < 7 || result.contains("#")){
				continue; //skip over non-result line
			}
			int scan = Integer.parseInt(tokens[1]);
			Scan s = parser.rap(scan);
			if(s.getHeader().getMsLevel() == 2){
				Spectrum spect = reader.getSpectrum(scan);
				spect.windowFilterPeaks2(15, 25);
				spect.mergePeaks(spect, 0.05);
				//spect.sqrtSpectrum();
				spect.peptide = tokens[4];
				int charge = Integer.parseInt(tokens[6]);
				TheoreticalSpectrum th = new TheoreticalSpectrum(spect.peptide+"."+charge);
				Spectrum projection = spect.project(th, 0.1);
				//Spectrum projection = spect;
				projection.computePeakRank();
				double score = scorer3.compare(th, projection);
				double[] stats = th.analyzeAnnotation(projection, spect.peptide+"."+charge, 0.03, false);
				
				System.out.print(result+"\t");
				System.out.print(score + "\t");
				for(int j = 0; j < stats.length; j++){
					System.out.print(stats[j] +"\t");
				}
				System.out.println();
				//System.out.println("mappped RT: " + spect.spectrumName + "\t" + spect.peptide + "\t" + SWATHUtils.getRT(s));
				//timeMap.put(SWATHUtils.getRT(s), spect);
			}
		}
	}
	
	/**
	 * This is an experimental method still under development.  In addition to spectral similarity measures
	 * this method also compute the match statistics usually used in PSM matches in database search (e.g. percent b/y, explained intensity etc...)
	 * SVM can then be used to learn a model that combined spectral-library searhc an database matches features to better filter out false positives
	 * @param minScan
	 * @param maxScan
	 */
	public static void testHybridLibrarySearch(int minScan, int maxScan){
		String ddaFile = "../mixture_linked/msdata/UPS_Ecoli//40fmol_UPS1_1ugEcoli_SWATH.mzXML";
		String ssmResult = "../mixture_linked/ACG_swathdevelopment_40fmol_UPS1_1ugEcoli_SWATH_vs_IDAlibrary_projcos_min0.6_top30_nodouble_pep.txt";
		MZXMLReader reader = new MZXMLReader(ddaFile);
		List<String> results = Utils.FileIOUtils.createListFromFile(ssmResult);
		
		//ConsensusSpectrumReader reader = new ConsensusSpectrumReader(queryFile);
		//reader.minInt=40;
		//reader.numNeighbors=3;
		Spectrum MS1 = null;
		SpectrumMap MSMap = null;
		SortedMap<Double, Spectrum> timeMap = new TreeMap<Double, Spectrum>();
		MSXMLParser parser = reader.getParser();
		String singlescorer = "../mixture_linked/yeast_NIST_lib_singlepep_win10_25.o";
		RankBaseScoreLearner pcomp = RankBaseScoreLearner.loadComparator(singlescorer);
		SimpleProbabilisticScorer scorer3 = new SimpleProbabilisticScorer(pcomp);
		scorer3.setMatchTolerance(0.05);
		for(int i = 0; i < results.size(); i++){
			String result = results.get(i);
			String[] tokens = result.split("\\t");
			if(tokens.length < 7 || result.contains("#")){
				continue; //skip over non-result line
			}
			int scan = Integer.parseInt(tokens[1]);
			Scan s = parser.rap(scan);
			if(s.getHeader().getMsLevel() == 2){
				Spectrum spect = reader.getSpectrum(scan);
				spect.windowFilterPeaks2(20, 25);
				spect.mergePeaks(spect, 0.05);
				//spect.sqrtSpectrum();
				spect.peptide = tokens[4];
				int charge = Integer.parseInt(tokens[6]);
				TheoreticalSpectrum th = new TheoreticalSpectrum(spect.peptide+"."+charge);
				Spectrum projection = spect.project(th, 0.1);
				//Spectrum projection = spect;
				projection.computePeakRank();
				double score = scorer3.compare(th, projection);
				double[] stats = th.analyzeAnnotation(projection, spect.peptide+"."+charge, 0.03, false);
				
				System.out.print(result+"\t");
				System.out.print(score + "\t");
				for(int j = 0; j < stats.length; j++){
					System.out.print(stats[j] +"\t");
				}
				System.out.println();
				//System.out.println("mappped RT: " + spect.spectrumName + "\t" + spect.peptide + "\t" + SWATHUtils.getRT(s));
				//timeMap.put(SWATHUtils.getRT(s), spect);
			}
		}
	}
	
	/**
	 * test incorporating the time-correlation similarity measure into the scoring function
	 * @param minScan
	 * @param maxScan
	 */
	public static void test2DSim(int minScan, int maxScan){
		String ddaFile = "../mixture_linked/msdata/UPS_Ecoli_Wiff/Duplicate_runs_201308/REP2/18488_REP3_40fmol_UPS1_1ug_Ecoli_NewStock2_SWATH_1.mzXML";
		String ssmResult = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPS_Ecoli_REP3withnewstock2_msplit_2_sortedscans.txt";
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_UPS12Ecoli_IDA_combined_RTlib_plusDecoy2.mgf";
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		MZXMLReader reader = new MZXMLReader(ddaFile);
		List<String> results = Utils.FileIOUtils.createListFromFile(ssmResult);
		ConsensusSpectrumReader reader2 = new ConsensusSpectrumReader(ddaFile);
		//reader.minInt=40;
		//reader.numNeighbors=3;
		Spectrum MS1 = null;
		SpectrumMap MSMap = null;
		SortedMap<Double, Spectrum> timeMap = new TreeMap<Double, Spectrum>();
		MSXMLParser parser = reader.getParser();
		for(int i = 0; i < results.size(); i++){
			String result = results.get(i);
			String[] tokens = result.split("\\t");
			if(tokens.length < 7 || result.contains("#")){
				continue; //skip over non-result line
			}
			if(tokens[4].contains("+0")){
				continue;
			}
			int scan = Integer.parseInt(tokens[1]);
			
			Scan s = parser.rap(scan);
			if(s.getHeader().getMsLevel() == 2){
				Spectrum query = reader.getSpectrum(scan);
				query.windowFilterPeaks2(15, 25);
				query.mergePeaks(query, 0.05);
				query.sqrtSpectrum();
				Spectrum libSpect = lib.getSpectra(tokens[4]+"."+tokens[6]).get(0);
				libSpect.mergePeaks(libSpect, 0.05);
				libSpect.filterPeaks(30);
				libSpect.sqrtSpectrum();	
				if(!tokens[4].equals("GAEPSGGAAR")){
					//continue;
				}
				if(Double.parseDouble(tokens[7]) < 0.75){
					//continue;
				}
				List<Spectrum> neighs = reader2.getNeighborScans(scan, 5, 5);
				double[] simTwoD = reader2.getProjectCosine(libSpect, neighs, 50);		
				double score = reader2.getProjectCosine(libSpect, query);
				double projCos = libSpect.projectedCosine(query, 0.05);
				libSpect.filterPeaks(10);
				double projCosTop = libSpect.projectedCosine(query, 0.05);
				System.out.println(tokens[0] + "\t" + tokens[1] + "\t" + tokens[4] + "\t" + tokens[8] + "\t" 
				+ tokens[7] + "\t" + tokens[11] +"\t"+ projCos + "\t" + projCosTop + "\t" + score 
				+ "\t" + simTwoD[0] + "\t" + simTwoD[1] + "\t" + simTwoD[2] 
			    + "\t" + simTwoD[3] + "\t" + simTwoD[4] + "\t" + simTwoD[5] + "\t" + simTwoD[6] + "\t" + simTwoD[7]
				+ "\t" + simTwoD[8] +"\t" + simTwoD[9] 
			    + "\t" + simTwoD[10] + "\t" + simTwoD[11]);
				//System.out.println(tokens[0] + "\t" + tokens[1] + "\t" + tokens[4] + "\t" + tokens[8] + "\t" + tokens[7] + "\t" + score + "\t" + projCos);
			}
		}
	}
	
	/**
	 * Given a list of candidate spectrum, return all SWATH spectrum that pass a minCos threshold
	 * @param specList
	 * @param s
	 * @param MSMap
	 * @param minCos
	 * @param parentMassTol
	 * @param fragMassTol
	 * @return
	 */
	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, SpectrumMap MSMap, double minCos, double parentMassTol, double fragMassTol){
		return bestPsimSpec(specList, s, MSMap, minCos, parentMassTol, fragMassTol, Double.MAX_VALUE);
	}
	
	/**
	 * Same as above but only consider candidates that have RT that is at most maxRTDiff away from query spectrum
	 * @param specList
	 * @param s
	 * @param MSMap
	 * @param minCos
	 * @param parentMassTol
	 * @param fragMassTol
	 * @param maxRTDiff
	 * @return
	 */
	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, SpectrumMap MSMap, double minCos, double parentMassTol, double fragMassTol, double maxRTDiff){
		return bestPsimSpec(specList, s, MSMap, minCos, parentMassTol, fragMassTol, maxRTDiff, Mass.DIFF_DA);
	}
	
	/**
	 * @param specList
	 * @param s
	 * @param MSMap
	 * @param minCos
	 * @param parentMassTol
	 * @param fragMassTol
	 * @param maxRTDiff
	 * @param mode - the fragment mass tolerance mode, can be either Da or PPM
	 * @return
	 */
	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, SpectrumMap MSMap, double minCos, double parentMassTol, double fragMassTol, double maxRTDiff, int mode){
		double left= s.parentMass-0.25*parentMassTol;
		double right=s.parentMass+0.9*parentMassTol;
		return bestPsimSpec(specList, s, MSMap, minCos, left, right, fragMassTol, maxRTDiff, mode);
	}
	
	/**
	 * Instead of specify an error tolerance, specify the actual parentmass interval [left, right] where to consider candidate spectrum
	 * @param specList
	 * @param s
	 * @param MSMap
	 * @param minCos
	 * @param left
	 * @param right
	 * @param fragMassTol
	 * @param maxRTDiff
	 * @param mode - fragment tolerance mode either in PPM or Da
	 * @return
	 */
	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, SpectrumMap MSMap, double minCos, double left, double right, double fragMassTol, double maxRTDiff, int mode){
		TreeMap<Double, Spectrum> bestList = new TreeMap();
		double bestScore = -1000000.0, currScore = 0.0;
		int count=0;
		Iterator<Spectrum> specIter = specList.iterator();
		List<Double> redundance = new ArrayList<Double>();
		int candCount=0;
		while(specIter.hasNext()){
			Spectrum s1 = specIter.next();
			if(s1.peptide.contains("0.9")){
				continue;
			}
			double offSetC13 = (Mass.C13-Mass.C12)/s1.charge;
//			if(Math.abs(s1.parentMass - s.parentMass) < parentMassTol 
//					|| Math.abs(s1.parentMass - offSetC13 - s.parentMass) < parentMassTol
//					|| Math.abs(s1.parentMass + offSetC13 - s.parentMass) < parentMassTol){
//			if(Math.abs(s1.parentMass - s.parentMass) < 3){
//			if(s1.parentMass > s.parentMass - 5 &&  
//					(s1.parentMass  - s.parentMass  < parentMassTol - 5)){
			double libRT = 0.0;
			if(s1.parentMass >= left &&  
						s1.parentMass  <= right  && 
						(Math.abs(s1.rt - s.rt ) < maxRTDiff)){

				//if(MSMap.checkPeak(s1.parentMass, 0.03)){
					//Spectrum s2 = new Spectrum(s1);
					//for(int i = 0; i < s1.getPeak().size(); i++){
					//	s2.getPeak().get(i).setIntensity(10);
					//}
					//s1.removePeaksInMass(s.parentMass-4, s.parentMass+21);
					currScore = s1.projectedCosine(s, fragMassTol, mode);
					//currScore = s2.projectedCosine(s, fragMassTol);
					//currScore = s1.projectedCosineWithSkip(s, fragMassTol, 1);
					
				//}
				candCount++;
			}else{
				currScore = 0;
			}
			if(currScore > minCos){
				Iterator<Double> it = bestList.descendingKeySet().iterator();
				boolean isFound = false;
				while(it.hasNext()){
					double currBest = it.next();
					Spectrum cand = bestList.get(currBest);
					//System.out.println("self-match " + cand.peptide + "\t" + s1.peptide + "\t" + 
					//+ cand.projectedCosine(s1, fragMassTol) + "\t" + s1.projectedCosine(cand, fragMassTol));
					//double shared = cand.sharePeaks(s1, fragMassTol);
					if((cand.projectedCosine(s1, fragMassTol, mode) > maxSelfCos
							|| s1.projectedCosine(cand, fragMassTol, mode) > maxSelfCos)){
						isFound = true;						
						if(Math.abs(currScore  - currBest) < 0.05){ //when indistinguishable, be conservative and don't call mod-peptide ID
							String seq1 = Utils.StringUtils.getStrippedSeq(s1.peptide);
							String seq2 = Utils.StringUtils.getStrippedSeq(cand.peptide);
							if(seq1.equals(seq2)){
								if(s1.peptide.contains("+")){
									break;
								}
								if(cand.peptide.contains("+")){
									redundance.add(currBest);
									isFound=false;
									break;
								}
							}
						}
						if(currScore > currBest){
							redundance.add(currBest);
							isFound = false;
						}
//						//System.out.println("find some dependency: " + s1.peptide + "\t" + cand.peptide);
						break;
					}
				}
				if(!isFound){
					bestList.put(currScore, s1);
				}
			}			
		}
		//System.out.println("Number of candidates: " + candCount);
		for(int i = 0; i < redundance.size(); i++){
			bestList.remove(redundance.get(i));
		}
		//System.out.println("done one");
		return bestList;
	}
	
	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, double minCos, double parentMassTol, double fragMassTol){
		TreeMap<Double, Spectrum> bestList = new TreeMap();
		TreeMap<Double, Spectrum> bestList2 = new TreeMap();
		double bestScore = -1000000.0, currScore = 0.0;
		int count=0;
		Iterator<Spectrum> specIter = specList.iterator();
		while(specIter.hasNext()){
			Spectrum s1 = specIter.next();
			double offSetC13 = (Mass.C13-Mass.C12)/s1.charge;
			if(Math.abs(s1.parentMass - s.parentMass) < parentMassTol 
					|| Math.abs(s1.parentMass - offSetC13 - s.parentMass) < parentMassTol
					|| Math.abs(s1.parentMass + offSetC13 - s.parentMass) < parentMassTol
			){
				//Spectrum s2 = new Spectrum(s1);
				//for(int i = 0; i < s1.getPeak().size(); i++){
				//	s2.getPeak().get(i).setIntensity(10);
				//}
				currScore = s1.projectedCosine(s, fragMassTol);
				//currScore = s2.projectedCosine(s, fragMassTol);
				//currScore = s1.projectedCosineWithSkip(s, fragMassTol, 1);
			}else{
				currScore = 0;
			}
			if(currScore > minCos){
				bestList.put(currScore, s1);
			}
						
		}
		//try interatively remove dependency
		Iterator<Double> it = bestList.descendingKeySet().iterator();
		Spectrum swath = s;
//		for(;it.hasNext();){
//			double score = it.next();
//			Spectrum cand = bestList.get(score);
//			double share = cand.sharePeaks(swath, fragMassTol);
//			double scoreNew = cand.projectedCosine(swath, fragMassTol);
//			//System.out.println("score difference: " + (scoreNew - score));
//			if(scoreNew > score){
//				System.out.println("improving");
//				bestList2.put(scoreNew, cand);
//			}
//			if(share > 10){
//				swath = s.removeSharePeaks(cand, fragMassTol);
//			}
//		}
		
		//System.out.println("size: " + bestList.size());
		bestList.putAll(bestList2);
		return bestList;
	}
	
	
	protected static TreeMap<Double,Spectrum> bestCondPsimSpec(TreeMap<Double, Spectrum> first, Spectrum s, double minCos, double fragMassTol){
		TreeMap<Double, Spectrum> bestList = new TreeMap();
//		double bestScore = -1000000.0, currScore = 0.0;
//		int count=0;
//		s.sqrtSpectrum();
//		Iterator<Spectrum> specIter = specList.iterator();
//		while(specIter.hasNext()){
//			Spectrum s1 = specIter.next();
//			double offSetC13 = (Mass.C13-Mass.C12)/s1.charge;
//			if(Math.abs(s1.parentMass - s.parentMass) < parentMassTol 
//					|| Math.abs(s1.parentMass - offSetC13 - s.parentMass) < parentMassTol
//					|| Math.abs(s1.parentMass + offSetC13 - s.parentMass) < parentMassTol
//					
//			){
//				currScore = s1.projectedCosine(s, fragMassTol);
//			}else{
//				currScore = 0;
//			}
//			if(currScore > minCos){
//				bestList.put(currScore, s1);
//			}
//		}
		return bestList;
	}
	
	

	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, int topN, double parentMassTol, double fragMassTol){
		return bestPsimSpec(specList.iterator(), s, topN, parentMassTol, fragMassTol);
	}
	
	/**
	 * Instead of specifiying a minimum cosine similarity for matches, this method specify the top N matches with highest similarity
	 * @param specIter
	 * @param s
	 * @param topN
	 * @param parentMassTol
	 * @param fragMassTol
	 * @return
	 */
	public static TreeMap<Double,Spectrum> bestPsimSpec(Iterator<Spectrum> specIter, Spectrum s, int topN, double parentMassTol, double fragMassTol){
		TreeMap<Double, Spectrum> bestList = new TreeMap();
		double bestScore = -1000.0, currScore = 0.0;
		int count=0;
		//s.sqrtSpectrum();
		count++;
		while(specIter.hasNext()){
			Spectrum s1 = specIter.next();
			if(Math.abs(s1.parentMass - s.parentMass) < parentMassTol){
				currScore = s1.projectedCosine(s, fragMassTol);
				count++;
				//currScore = s1.cosineApprox(s, fragMassTol);
				//currScore = s1.cosine(s, fragMassTol);
			}else{
				currScore = -100000;
			}
			if(currScore > bestScore){
				bestList.put(currScore, s1);
				if(bestList.size() > topN){
						bestScore = bestList.firstKey().doubleValue();
						bestList.remove(bestScore);
						bestScore = bestList.firstKey().doubleValue();
					}
				}
			}
		//System.out.println("number of candidates: " + count);
		//System.out.println("bestList: " + bestList.size());
		return bestList;
	}
	
	

	
	protected static void getSWATHSpectrum(int index, String spectFile){
		MZXMLReader reader = new MZXMLReader(spectFile);
		Spectrum s = reader.getSpectrum(index);
		System.out.println(s.toString());
	}
	
	
	/**
	 * The SWATH window widths can be specified as a input file, whcih specify
	 * the particular acquisition methods.  This is helpful when using multiple
	 * different window widths during a acquisiton cycle. Parse the window information from 
	 * input file
	 * @param file
	 * @return
	 */
	public List<double[]> parseSWATHWindows(String file){
		List<String> lines = Utils.FileIOUtils.createListFromFile(file);
		List<double[]> windows = new ArrayList<double[]>(lines.size());
		System.out.println("parsing swath window file");
		int cycle = 0;
		for(int i = 0; i < lines.size(); i++){
			System.out.println(lines.get(i));
			String[] tokens = lines.get(i).split("\\s+");
			if(lines.get(i).startsWith("#")){
				continue;
			}
			double[] window = new double[2];
			try{
				if(tokens.length >= 3){
					if(tokens[0].equals("MS2")){
						window[0] = Double.parseDouble(tokens[1]);
						window[1] = Double.parseDouble(tokens[2]);
						windows.add(window);
					}
					cycle++;
				}
			}catch(NumberFormatException e){
				System.err.println("Invalid format in the swath window file");
				e.printStackTrace();
				System.exit(1);
			}
		}
		if(windows.size() > 0){
			System.out.println("Parsed " + windows.size() + " isolation windows from file: " + file);
			this.swathwindows = windows;
			this.SWATHCycle = cycle;
			System.out.println("Cycle: " + cycle);
		}
		return windows;
	}
	
	protected static void testCrossLibrarySimilarity(){
		String libFile1 = "../mixture_linked/gringar_swath_SLRIACG01_msgfdb_ids.mgf";
		String libFile2 = "../mixture_linked/gringar_SLRIACG01_cdk47ct_msgfdbids.mgf";
		SpectrumLib lib1 = new SpectrumLib(libFile1, "MGF");
		SpectrumLib lib2 = new SpectrumLib(libFile2, "MGF");
		List<Spectrum> specList = lib1.getAllSpectrums();
		for(int i = 0; i < specList.size(); i++){
			Spectrum qtof = specList.get(i);
			List<Spectrum> ionTraps = lib2.getSpectra(qtof.peptide+"."+qtof.charge);
			//System.out.println("current peptide: " + qtof.peptide);
			if(ionTraps != null){
				for(int j = 0; j < ionTraps.size(); j++){
					Spectrum IT = ionTraps.get(j);
					IT.sqrtSpectrum();
					qtof.sqrtSpectrum();
					IT.shiftSpectrumPPM(100);
					//qtof.shiftSpectrumPPM(10);	
					IT.windowFilterPeaks2(6, 25);
					qtof.windowFilterPeaks2(6, 25);
					//IT.filterPeaks(20);
					//qtof.filterPeaks(20);
					if(!(IT.scanNumber == 4797 &&  qtof.scanNumber == 4513)){
						//continue;
					}
					//IT.filterPeaks(20);
					//qtof.filterPeaks(20);
					double similarity = IT.cosine(qtof, 0.1);
					double count = IT.sharePeaks(qtof, 0.1);
					double psim1 = IT.projectedCosine(qtof, 0.1);
					double psim2	 = qtof.projectedCosine(IT, 0.1);
					System.out.println("Similarity: " + qtof.spectrumName + "\t" + qtof.peptide+ "\t" + qtof.charge 
							+ " and " + IT.spectrumName + "\t" + IT.peptide + "\t" + IT.peptide +"\t" + similarity 
							+ "\t" + psim1 + "\t" + psim2 +"\t" + count);
				}
			}
		}
	}
	

	protected static void simpleCosTest(){
		String queryFile = "../mixture_linked/testlib.mgf";
		String queryFile2 = "../mixture_linked/msdata/gringar/cdk47ct_swath-cdk47ct.mzXML";
		//MZXMLReader reader = new MZXMLReader(queryFile);
		SpectrumLib lib = new SpectrumLib(queryFile, "MGF");
		MZXMLReader reader2 = new MZXMLReader(queryFile2);
		//Spectrum s1 = reader.getSpectrum(1631); 
		System.out.println("lib size: " + lib.getSpectrumList().size());
		Spectrum s1 = lib.getSpectrumList().get(0);
		Spectrum s2 = reader2.getSpectrum(5000);
		s1.sqrtSpectrum();
		s2.sqrtSpectrum();
		s1.windowFilterPeaks(6, 25);
		s2.windowFilterPeaks2(20, 25);
		s2.shiftSpectrumPPM(100);
		System.out.println("matching: " + s1.spectrumName + "\t" + s1.parentMass + "\t" + s1.charge 
				+ " to: " + s2.spectrumName + "\t" + s2.parentMass + "\t" + s2.charge);
		System.out.println(" cosine: " + s1.projectedCosine(s2, 0.05));
	}
	
	
	
	
	public static void testSearch(String[] args){
		String queryFile = "../mixture_linked//msdata/UPS_Ecoli_Wiff/Duplicate_runs_201308/REP2/18482_REP3_500ng_Ecoli_NewStock2_SWATH_2.mzXML";
		String outFile = "../mixture_linked/test_EmilyToni_withNISTLib.txt";
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_UPSEcoli_REP234_IDA_plusDecoy2.mgf";
//		String modFile = "../mixture_linked/Mod_MSPLIT-DIA.txt";
//		String swathWindowsFile = "../mixture_linked\\test_swath_windows.txt";
		args = new String[]{"25", "50", "0", queryFile, libraryFile, outFile};
		CommandLineParser cmdParser = new CommandLineParser(args);
		SWATHMSPLITSearch search = new SWATHMSPLITSearch();
		search.parent = cmdParser.getDouble(0);
		search.fragment = cmdParser.getDouble(1);
		search.queryFile = cmdParser.getString(3);
		search.SWATHCycle = cmdParser.getInteger(2);
		if(args.length > 6){
			String ptmFile = cmdParser.getString(6);
			List<PTM[]> ptmList=null;
			if(ptmFile.endsWith(".xml")){
				ptmList = PTM.parsePTMFromXML(ptmFile);
			}else{
				ptmList = PTM.parsePTMs(ptmFile);
			}
			search.ptmList = ptmList;
		}
		if(args.length > 7){
			String windowsFile = cmdParser.getString(7);
			search.parseSWATHWindows(windowsFile);
		}
		if(search.SWATHCycle <= 0){
			search.SWATHCycle = search.getSWATHCycle();//cmdParser.getInteger(2);
		}
		search.libraryFile = cmdParser.getString(4);
		String out = cmdParser.getString(5);
		search.outFile = Utils.FileIOUtils.generateOutFile(out, search.queryFile, "_msplitout.txt");
		search.startMSPLITSearch(0, 1000000000);
	}
	
	public static void main(String[] args){
		testSearch(args);
		//targetedIdentification();
		//testMSPLITSearch(10, 1000000, Double.parseDouble(args[0]), Double.parseDouble(args[1]), args[2], args[3], args[4]);
		//testReverseMSPLITSearch(100, 100000);
		//testHybridLibrarySearch(10,100000);

		//testCrossLibrarySimilarity();
		
	}

}
