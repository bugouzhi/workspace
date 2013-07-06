package org.Spectrums;

import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Initial try on searching SWATH data with M-SPLIT
 * @author Jian Wang	
 *
 */
public class SWATHMSPLITSearch {
	public static void testMSPLITSearch(int minScan, int maxScan){
		String queryFile = "../mixture_linked/msdata/gringar/swath_development/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		//String libraryFile = "../mixture_linked/NIST_human_IT_2010-01-14_plusDecoy.mgf";
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
		Iterator<Spectrum> reader = new MZXMLReader(queryFile);
		//ConsensusSpectrumReader reader = new ConsensusSpectrumReader(queryFile);
		//reader.minInt=40;
		//reader.numNeighbors=3;
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		lib.windowFilterPeaks(6, 25);
		lib.mergeSpectrum(0.075);
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
		System.out.println("#File\tScan#\tMz\tz\tPeptide\tMz\tz\tcosine\tName\t#Peak(Query)\t#Peaks(match)\t#shared\tfraction-matched\trelative-alpha\tIonCount");
		while(iter.hasNext()){
			Spectrum s = iter.next();
			minScan = 1000;
			//maxScan = 41;
			if(s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}else{
				if(go){
					start = (new GregorianCalendar()).getTimeInMillis();
					go = false;
				}				
			}
			//s.filterPeaksByIntensity(100);
			//System.out.println(s);
			//s.toRelIntensity();
			//s.computePeaksZScore(0.9);
			//System.out.println(s);
			s.filterPeaksByIntensity(50);
			//s.filterPeaksByRankScore(6);
			s.filterPeaks(500);
			//s.windowFilterPeaks2(15, 25);
			s.mergePeaks(s, 0.075);
			//System.out.println(s);
			//System.out.println("Spectrum: " + s.scanNumber + "\t" + s.parentMass +"\t" + s.charge + "\t#peaks: " + s.getPeaks().size());
			//TreeMap<Double,Spectrum> bestCands = bestPsimSpec(lib.getAllSpectrums(), s, 1, 2, 0.05);
			//s.shiftSpectrumPPM(100);
			TreeMap<Double,Spectrum> bestCands = bestPsimSpec(lib.getAllSpectrums(), s, 0.6, 25, 0.05);
			//TreeMap<Double,Spectrum> bestCands = bestPsimSpec(reader, s, 10, 25, 0.03);
//			System.out.println("best cand size: " + bestCands.size());
			//System.out.println("best score: " + bestCands.lastKey());
			Iterator it = bestCands.descendingKeySet().iterator();
			double maxInt = 0.0; //matches with maximum abundance
			while(it.hasNext()){
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				double projectInt = cand.projectedPeakIntensity(s, 0.05);
				maxInt = projectInt > maxInt ? projectInt : maxInt;
			}
			
			it = bestCands.descendingKeySet().iterator();
			while(it.hasNext()){	
				Double psim = (Double)it.next();
				Spectrum cand = bestCands.get(psim);
				if(cand.protein == null){
					cand.protein = cand.spectrumName;
				}
				double sharePeaks = cand.sharePeaks(s, 0.1);
				double projectInt = cand.projectedPeakIntensity(s, 0.05);
				double projectInt2 = s.projectedPeakIntensity(cand, 0.05);
				if(sharePeaks > 0){
					System.out.println(queryFile + "\t" +  s.scanNumber + "\t" + s.parentMass +"\t" + s.charge +"\t" 
							+ cand.peptide + "\t"  + cand.parentMass + "\t" + cand.charge + "\t" + psim 
							+ "\t" + cand.protein + "\t" + s.getPeaks().size() + "\t" + cand.getPeak().size() 
							+"\t" + sharePeaks + "\t" + projectInt + "\t" + projectInt/maxInt + "\t" + s.upperBound +"\t" );
					//System.out.println(s);
					//System.out.println();
					//System.out.println(cand);
				}
			}
			counter++;
			if(counter == 5){
				//break;
			}
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	
	public static void testReverseMSPLITSearch(int minScan, int maxScan){
		String queryFile = "../mixture_linked/msdata/gringar/swath_development/14341_UPS1-400fm_IDA_5600 .mzXML";
		//String libraryFile = "../mixture_linked/spectral_library/human_qtof_plusDecoy.sptxt";
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
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
			s1.windowFilterPeaks(6, 25);
			s1.mergePeaks(s1, 0.075);
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
			s.filterFlatPeaks(0.05);
			s.windowFilterPeaks2(15, 25);
			s.mergePeaks(s, 0.075);
			s.sqrtSpectrum();
			System.out.println("Spectrum: " + s.scanNumber + "\t" + s.parentMass +"\t" + s.charge + "\t#peaks: " + s.getPeaks().size());
			
			for(int i = 0; i < specList.size(); i++){
				Spectrum s1 = specList.get(i);
				if(Math.abs(s1.parentMass - s.parentMass) < parentMassTol){
					double currScore = s1.projectedCosine(s, 0.05);
					if(scoreMap.get(s1) < currScore){
						scoreMap.put(s1, currScore);
						bestMap.put(s1, s);
					}
					
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
				System.out.println(s1.spectrumName + "\t" + s1.peptide + "\t" + s1.parentMass + "\t" + s1.charge 
					+ "\t" +  best.scanNumber + "\t" + score + "\t" + best.getPeak().size() + "\t" + s1.getPeak().size() + "\t" + sharePeaks);
			}
		}

		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	
	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, double minCos, double parentMassTol, double fragMassTol){
		TreeMap<Double, Spectrum> bestList = new TreeMap();
		double bestScore = -1000000.0, currScore = 0.0;
		int count=0;
		s.sqrtSpectrum();
		Iterator<Spectrum> specIter = specList.iterator();
		while(specIter.hasNext()){
			Spectrum s1 = specIter.next();
			if(Math.abs(s1.parentMass - s.parentMass) < parentMassTol){
				currScore = s1.projectedCosine(s, fragMassTol);
			}else{
				currScore = 0;
			}
			if(currScore > minCos){
				bestList.put(currScore, s1);
			}
		}
		return bestList;
	}
	
	public static TreeMap<Double,Spectrum> bestPsimSpec(List<Spectrum> specList, Spectrum s, int topN, double parentMassTol, double fragMassTol){
		return bestPsimSpec(specList.iterator(), s, topN, parentMassTol, fragMassTol);
	}
	
	public static TreeMap<Double,Spectrum> bestPsimSpec(Iterator<Spectrum> specIter, Spectrum s, int topN, double parentMassTol, double fragMassTol){
		TreeMap<Double, Spectrum> bestList = new TreeMap();
		double bestScore = -1000.0, currScore = 0.0;
		int count=0;
		s.sqrtSpectrum();
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
	
	

	
	public static void getSWATHSpectrum(int index, String spectFile){
		MZXMLReader reader = new MZXMLReader(spectFile);
		Spectrum s = reader.getSpectrum(index);
		System.out.println(s.toString());
	}
	
	public static void testCrossLibrarySimilarity(){
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
	

	public static void simpleCosTest(){
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
	
	public static void main(String[] args){
		testMSPLITSearch(10, 1000000);
		//testReverseMSPLITSearch(100, 100000);
		//getSWATHSpectrum(5000, "../mixture_linked/msdata/gringar/swath_development/14345_UPS2-4pm_IDA_5600.mzXML");
		//simpleCosTest();
		//testCrossLibrarySimilarity();
	}

}
