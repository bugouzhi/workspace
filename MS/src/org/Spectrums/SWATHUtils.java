package org.Spectrums;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.systemsbiology.jrap.stax.Scan;

/**
 * Some utility function for SWATH data
 * @author Jian Wang
 *
 */
public class SWATHUtils {
	public static final int DALTON = 1;
	public static final int PPM = 2;
	
	
	public static boolean checkPrecursor(double ms2, double precursor, double windowWidth){
		 double diff = precursor - ms2 + 4.0;
		 return diff > 0 && diff < windowWidth;
	}
	
	public static boolean checkMass(double ms2, double precursor, double tolerance, int mode){
		return SpectrumUtil.checkMass(ms2, precursor, tolerance, mode);
	}
	
	
	public static int[] peakIntDistr(Spectrum s, double min, double max, double binWidth){
		List<Peak> sorted = new ArrayList();
		sorted.addAll(s.getPeak());
		Collections.sort(sorted, PeakIntensityComparator.comparator);
		int bins = (int)Math.ceil((max - min)/binWidth);
		int[] counts = new int[bins];
		for(int i = 0; i < sorted.size(); i++){
			
		}
		return null;
	}
	
	public static double getRT(Scan s){
		String rt = s.getHeader().getRetentionTime();
		return Double.parseDouble(rt.substring(2, rt.length()-1));
	}
	
	/**
	 * get the corresponding swath ms1 scan
	 * @param s
	 * @return
	 */
	public static int getSWATHMS1Scan(Spectrum s){
		return (int)(s.scanNumber / 35.0)*35+1;
	}
	
	//convert a spectral library to the peakview ion library text format
	public static void getIonLibrary(String spectrumFile){
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		List<Spectrum> specList = lib.getAllSpectrums();
		getIonLibrary(specList);
	}
	
	public static void getIonLibrary(List<Spectrum> specList){
		int counter = 0;
		int ProteinCounter = 1; //not sure if this is protein group number but in ionlibrary seems unique to each proteins
		String header = "Q1\tQ3\tRT_detected\tprotein_name\tisotype\trelative_intensity\tstripped_sequence\tmodification_sequence\tprec_z\tfrg_type\tfrg_z\tfrg_nr\tiRT\tuniprot_id\tscore\tdecoy\tprec_y\tconfidence\tshared\tN\tmods\tnterm\tcterm";
		System.out.println(header);
		Map<String, Integer> proteinMap = new HashMap<String, Integer>();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.peptide.contains("+") || s.spectrumName.contains("DECOY")){
				continue;
			}
			SimpleMatchingGraph g = SpectrumUtil.getMatchGraph(s, 0.05);
			TreeSet sortedSet = new TreeSet(PeakIntensityComparator.comparator);
			sortedSet.addAll(g.vertexSet(SimpleMatchingGraph.Observed));
			int matchCount = 0;
			StringBuffer outbuffer = new StringBuffer();
			//for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			for(Iterator<LabelledPeak> it = sortedSet.descendingIterator(); it.hasNext();){
				Peak p = it.next();
				List neighs = g.getNeighbors(p);
				LabelledPeak closestP = null;
				double smallestErr = 10000;
				for(int n = 0; n < neighs.size(); n++){
					LabelledPeak currPeak = (LabelledPeak)neighs.get(n);
					double currDiff = Math.abs(currPeak.getMass()-p.getMass());
					closestP = currDiff < smallestErr ? currPeak : closestP;
					smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
				}
				int protInd = 1;
				
				String[] tokens = s.spectrumName.split("\\s+");
				if(s.rt == 0){
					s.rt = Double.parseDouble(tokens[5].substring(2, tokens[5].length()-1));
				}
				if(s.protein == null){
					s.protein = tokens[7];	
				}
				if(proteinMap.containsKey(s.protein)){
					protInd = proteinMap.get(s.protein);
				}else{
					protInd = ProteinCounter++;
					proteinMap.put(s.protein, protInd);
				}
				double Int = 10000;
				if(closestP != null 
						//&& (closestP.getType().equals("b") || closestP.getType().equals("y"))
						&&		closestP.getPos() > 1){
					matchCount++;
					outbuffer.append(s.parentMass + "\t");
					outbuffer.append(p.getMass() + "\t");
					outbuffer.append(s.rt/60 + "\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(" " + "\t");
					outbuffer.append(p.getIntensity() + "\t");
					outbuffer.append(s.peptide.replaceAll("[0-9\\.\\+\\-]", "") +"\t");
					outbuffer.append(s.peptide.replaceAll("[0-9\\.\\+\\-]", "") +"\t");
					//System.out.print(s.peptide + "\t");
					outbuffer.append(s.charge  + "\t");
					outbuffer.append(closestP.getType() + "\t");
					outbuffer.append(closestP.getCharge() + "\t");
					outbuffer.append(closestP.getPos() + "\t");
					outbuffer.append(s.rt/60 + "\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(-1 + "\t");
					outbuffer.append("False" + "\t");
					outbuffer.append(Int + "\t");
					outbuffer.append("1" + "\t");
					outbuffer.append("False" + "\t");
					outbuffer.append(protInd + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("\n");
				}
				
				if(matchCount > 500 ){
					break;
				}

			}
			if(matchCount >= 5){
				System.out.print(outbuffer);
			}else{
				//System.out.println("warning less than 10 peaks in library spectrum");
			}
			counter++;
		}
	}
	
	/**
	 * Generate a ion library for targeted extraction based on identification from multiplexed 
	 * search tools such as M-SPLIT, the identified peptides are then used to do targeted extraction from swath
	 * @param libFile
	 * @param SpectrumFile
	 * @param resultsFile
	 */
	public static void getIonLibrary(String libFile, String spectrumFile, String resultsFile, String ssmResult){
		SpectrumLib lib = new SpectrumLib(libFile, "MGF");
		ProteinIDExtractor protID = new ProteinIDExtractor(lib.getAllSpectrums(), "../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta");
		Map<String, List<String>> pMap = protID.peptideMap;
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		SearchResult searchResult = new SearchResult(ssmResult, "", lib, reader);
		searchResult.parseResult(0, 1, 4, 6, 8, 31);
		List<String> results = Utils.FileIOUtils.createListFromFile(resultsFile);
		List<Spectrum> specList = new ArrayList<Spectrum>();
		int counter = 0;
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String[] result = it.next().split("\\t");
			if(result.length < 10){
				continue;
			}
			String pepID = result[4] + "." + result[6];
			if(lib.getSpectra(pepID) != null){
				Spectrum libSpect = lib.getSpectra(pepID).get(0);
				Spectrum swathSpect = reader.getSpectrum(Integer.parseInt(result[1]));
				//MixSSM mSSM = searchResult.getMatches(swathSpect.scanNumber);
				//if(mSSM == null){
				//	continue;
				//}
				//Spectrum unique = mSSM.getUniqueSpect(libSpect);
				//System.out.println("Spectrum size: " + libSpect.getPeak().size() + "\t" + unique.getPeak().size());
				libSpect.rt = swathSpect.rt;
				libSpect.scaleSpectrum(1/libSpect.sumMagnitude());
				//libSpect.setPeaks(unique.getPeak());
				//String[] tokens = libSpect.spectrumName.split("\\s+");
				//libSpect.rt = Double.parseDouble(tokens[5].substring(2, tokens[5].length()-1));
				
				if(pMap.containsKey(libSpect.peptide)){
					List<String> prots = pMap.get(libSpect.peptide);
					libSpect.score = prots.size();
					libSpect.protein = prots.get(0);
					specList.add(libSpect);
					counter++;
				}
				if(counter % 1000 == 0){
					System.out.println("got spec: " + counter);
				}
			}
		}
		System.out.println("Getting library for peptides: " + specList.size());
		getIonLibrary(specList);
	}
	
	/**
	 * Generate ion library base on fragments ions from swath directly
	 * (as oppose to those from reference library)
	 * @param libFile
	 * @param spectrumFile
	 * @param resultsFile
	 */
	public static void getIonLibraryFromSWATH(String libFile, String spectrumFile, String resultsFile){
		SpectrumLib lib = new SpectrumLib(libFile, "MGF");
		ProteinIDExtractor protID = new ProteinIDExtractor(lib.getAllSpectrums(), "../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta");
		Map<String, List<String>> pMap = protID.peptideMap;
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		List<String> results = Utils.FileIOUtils.createListFromFile(resultsFile);
		List<Spectrum> specList = new ArrayList<Spectrum>();
		int counter = 0;
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String[] result = it.next().split("\\t");
			String ID = result[4]+"."+result[6];
			if(lib.getSpectra(ID) != null){
				Spectrum libSpect = lib.getSpectra(ID).get(0);
				Spectrum swathSpect = reader.getSpectrum(Integer.parseInt(result[1]));
				//swathSpect.windowFilterPeaks2(15, 25);
				//swathSpect.mergePeaks(swathSpect, 0.05);
				Spectrum projection = swathSpect.project(libSpect, 0.05);
				projection.mergePeaks(projection, 0.05);
				//projection = projection.project(libSpect, 0.05);
				projection.rt = swathSpect.rt;
				projection.parentMass = libSpect.parentMass;
				projection.charge = libSpect.charge;
				projection.peptide = libSpect.peptide;
				//System.out.println(libSpect.peptide);
				if(pMap.containsKey(libSpect.peptide)){
					List<String> prots = pMap.get(libSpect.peptide);
					projection.score = prots.size();
					projection.protein = prots.get(0);
					specList.add(projection);
					counter++;
				}
				if(counter % 1000 == 0){
					System.out.println("got spec: " + counter);
				}
			}
		}
		//System.out.println(specList.size());
		getIonLibrary(specList);
	}
	
	public static void testGenerateIonLibrary(){
		String specLibFile = "../mixture_linked/APSWATH_PPP2RA_MSGFDB_IDs_1pepFDR_plusDecoy2.mgf";
		String spectrumFile = "J://workspace/mixture_linked/msdata/gringar/APSWATH/PPP2R1A/SWATH/19388_ZYW_emptyvector_PPP2R1A_rep2_SWATH_TOF1_chp_P107.mzXML";
		String resultFile = "../mixture_linked/t0";
		String ssmResultFile = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPSEcoli_libREp234_msplit_7_pep.txt";
		//getIonLibrary(specLibFile);
		getIonLibrary(specLibFile, spectrumFile, resultFile, ssmResultFile);
		//getIonLibraryFromSWATH(specLibFile, spectrumFile, resultFile);
	}
	public static void main(String[] args){
		testGenerateIonLibrary();
	}
}

