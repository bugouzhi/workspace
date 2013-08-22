package org.Spectrums;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A spectrum lib, indexed by peaks in the spectrum rather than
 * peptides
 * @author Jian Wang
 *
 */
public class LookUpSpectrumLib{
	private Map<Double, List<Spectrum>> table;
	List<Spectrum> specList;
	private double massTolerance;
	public LookUpSpectrumLib(List<Spectrum> specList){
		this.specList = specList;
		table = new HashMap<Double, List<Spectrum>>();
		System.out.println("indexing " + specList.size() + " spectra");
		indexSpectrumByPeak();
		//System.out.println("after indexing we have " + table.keySet().size() + " entries");
	}
	
	private void indexSpectrumByPeak(){
		for(Iterator<Spectrum> it = specList.iterator(); it.hasNext();){
			Spectrum s = it.next();
			List<Peak> peaks = s.getPeaks();
			//System.out.println("spectra has " + peaks.size() + "peaks");
			for(int i = 0; i < peaks.size(); i++){
				LabelledPeak lp = (LabelledPeak)peaks.get(i);
				if(lp.getCharge() > 1){
					continue;
				}
				Double key =new Double(Math.round(lp.getMass()));
				if(table.containsKey(key)){
					table.get(key).add(s);
				}else{	
					List<Spectrum> l = new ArrayList();
					l.add(s);
					table.put(key, l);
				}
			}
		}
	}
	
	public Set<Spectrum> getSpectrumByPeak(Peak p){
		Set<Spectrum> combine = new HashSet<Spectrum>();
		return getSpectrumByPeak(p, combine);
	}
	
	public Set<Spectrum> getSpectrumByPeak(Peak p, Set<Spectrum> resultSet){
		Double key = new Double(Math.round(p.getMass()));
		//account for mass error in the real spectrum
		Double key2 = new Double(Math.round(p.getMass() - massTolerance));
		Double key3 = new Double(Math.round(p.getMass() + massTolerance));
		if(table.containsKey(key)){
			resultSet.addAll(table.get(key));
		}
		if(!key.equals(key2) && table.containsKey(key2)){
			resultSet.addAll(table.get(key2));
		}
		if(!key.equals(key3) && !key2.equals(key3)&& table.containsKey(key3)){
			resultSet.addAll(table.get(key3));
		}
		return resultSet;
	}
	
	public Map<Spectrum, Integer> getSpectrumByPeak(Peak p, Map<Spectrum, Integer> resultMap){
		Double key = new Double(Math.round(p.getMass()));
		//account for mass error in the real spectrum
		Double key2 = new Double(Math.round(p.getMass() - massTolerance));
		Double key3 = new Double(Math.round(p.getMass() + massTolerance));
		if(table.containsKey(key)){
			incrementSpectrumCount(resultMap, table.get(key));
		}
		if(!key.equals(key2) && table.containsKey(key2)){
			incrementSpectrumCount(resultMap, table.get(key2));
		}
		if(!key.equals(key3) && !key2.equals(key3)&& table.containsKey(key3)){
			incrementSpectrumCount(resultMap, table.get(key3));
		}
		return resultMap;
	}
	
	private void incrementSpectrumCount(Map<Spectrum, Integer> m, List<Spectrum> s){
		for(Iterator<Spectrum> it = s.iterator(); it.hasNext();){
			Spectrum curr = it.next();
			if(m.containsKey(curr)){
				Integer count = m.get(curr);
				m.put(curr, count+1);
			}else{
				m.put(curr, new Integer(1));
			}
		}
	}
	
	public List<Spectrum> getSpectrumByPeaks(List<Peak> peakList){
		Set<Spectrum> combined = new HashSet<Spectrum>();
		List<Spectrum> l = new ArrayList<Spectrum>();
		for(int i = 0; i < peakList.size(); i++){
			getSpectrumByPeak(peakList.get(i), combined);
		}
		l.addAll(combined);
		return l;
	}
	
	public List<Spectrum> getSpectrumByPeaks(List<Peak> peakList, int minMatch){
		List<Spectrum> l = new ArrayList<Spectrum>();
		Map<Spectrum, Integer> hitsCount = new HashMap<Spectrum, Integer>();
		for(int i = 0; i < peakList.size(); i++){
			getSpectrumByPeak(peakList.get(i), hitsCount);
		}
		for(Iterator<Spectrum> it = hitsCount.keySet().iterator(); it.hasNext();){
			Spectrum s = it.next();
			if(hitsCount.get(s).intValue() > minMatch){
				l.add(s);
			}
		}
		return l;
	}
	
	public List<Spectrum> getSpectrumByPeaksStratified(List<Peak> peakList, int[] intervals, int[] minMatch){
		List<Spectrum> l = new ArrayList<Spectrum>();
		Set<Spectrum> specSet = new HashSet<Spectrum>();
		specSet.addAll(this.specList);
		int i = minMatch.length-1;
		do{
			Set<Spectrum> specSet2 = new HashSet<Spectrum>();
			Map<Spectrum, Integer> hitsCount = new HashMap<Spectrum, Integer>();
			Map<Spectrum, Integer> hitsCount2 = new HashMap<Spectrum, Integer>();
			for(int j = intervals[minMatch.length]; j >= intervals[i] && j > 0; j--){
				getSpectrumByPeak(peakList.get(j), hitsCount);
			}
			for(Iterator<Spectrum> it = hitsCount.keySet().iterator(); it.hasNext();){
				Spectrum s = it.next();
				if(hitsCount.get(s).intValue() > minMatch[i]){
					specSet2.add(s);
				}
			}
			getSetInterSect(specSet2, specSet);
			System.out.println("candidates so-far: " + specSet2.size());
			specSet = specSet2;
		}while(--i >=0);
		l.addAll(specSet);
		return l;
	}
	
	public void getSetInterSect(Set<?> set1, Set<?> set2){
		List toBeRemove = new ArrayList();
		for(Iterator it = set1.iterator(); it.hasNext();){
			Object o = it.next(); 
			if(!set2.contains(o)){
				toBeRemove.add(o);
			}
		}
		set1.removeAll(toBeRemove);
	}
	
	public static boolean checkPassFilter(String peptide, List<Spectrum> filtered){
		//System.out.println("checking peptide: " + peptide);
		for(int i = 0; i < filtered.size(); i++){
			Spectrum s = filtered.get(i);
			//System.out.println("peptide is: " + s.peptide);
			if(s.peptide.equals(peptide)){
				return true;
			}
		}
		return false;
	}
	
	public static void testLookUpSpectrum(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		String file = "..\\mixture_linked\\Ecoli_allpeptides.txt";
		lib1.removeModSpectra();
		List<Spectrum> specList = lib1.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < 1000; i++){
			Spectrum s = lib1.getRandomSpectrum();
			s.windowFilterPeaks(5, 25);
			s.computePeakRank();
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 9, false);
			System.out.println("After filter one we have: " + lib.getAllSpectrums().size() + " candidates");
			LookUpSpectrumLib filteredLib = new LookUpSpectrumLib(lib.getSpectrumList());
			List<Peak> queryPeaks = s.getTopPeaks(25);
			//List<Peak> queryPeaks = s.topWindowPeaks(6, 200);
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<Spectrum> candidates = filteredLib.getSpectrumByPeaks(queryPeaks, 4);
			boolean passedFilter = checkPassFilter(s.peptide, candidates);
			System.out.println("After filter two we have: " + candidates.size() + " candidates ");
			System.out.println("After filter correct peptide is retained?: " + passedFilter);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void testLookUpMixtureSpectrum(){
		String spectrumFile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib.txt";
		String mixturefile = "..\\mixture_linked\\exp_sim_mixtures.id"; 
//		SpectrumLib mixlib = lib2.createMix(mixturefile, 110, 1, 0.2, 0.0000000, 1.0, 3, false);
		SpectrumLib mixlib = lib2.createRandomMix(100, 1, 0.000003, 0.0000000, 1.0, 3, false);
		List<Spectrum>  specList = mixlib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		factory.setMinCharge(2);
		factory.setMaxCharge(3);
		System.out.println("Total memory used to index library: " + (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()));
		long start = (new GregorianCalendar()).getTimeInMillis();
		System.out.println("mixlib has size: " + specList.size());
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			LookUpSpectrumLib filteredLib = new LookUpSpectrumLib(lib.getSpectrumList());
			List<Peak> queryPeaks = s.getTopPeaks(26);
			System.out.println("After filter one we have: " + lib.getAllSpectrums().size() + " candidates");
			System.out.println("Query peaks has: " + queryPeaks.size());
//			List<Spectrum> candidates = filteredLib.getSpectrumByPeaksStratified(queryPeaks, new int[]{0, 5, 10, 15,20,25}, new int[]{4,3,2,1,0});
			List<Spectrum> candidates = filteredLib.getSpectrumByPeaks(queryPeaks, 2);
			String[] peps = s.peptide.split(" & ");
			boolean passedFilter = checkPassFilter(peps[0], candidates);
			boolean passedFilter2 = checkPassFilter(peps[1], candidates);
			System.out.println("After filter two we have: " + candidates.size() + " candidates ");
			System.out.println("After filter correct peptide is retained?: " + passedFilter + "\t" + passedFilter2);			
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	//we first find lookup score find first spectrum, then remove, lookup and find 2nd spectrum
	public static void testIterativeLookUpMixtureSpectrum(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";		
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib.txt";
		String mixturefile = "..\\mixture_linked\\exp_sim_mixtures.id"; 
		lib1.removeModSpectra();
		lib1.computeRank();
		SpectrumLib mixlib = lib2.createRandomMix(100, 1, 0.3, 0.0000000, 1.0, 3, false);
		List<Spectrum>  specList = mixlib.getSpectrumList();
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		factory.setMinCharge(2);
		factory.setMaxCharge(3);
		System.out.println("Total memory used to index library: " + (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()));
		long start = (new GregorianCalendar()).getTimeInMillis();
		System.out.println("mixlib has size: " + specList.size());
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			LookUpSpectrumLib filteredLib = new LookUpSpectrumLib(lib.getSpectrumList());
			List<Peak> queryPeaks = s.getTopPeaks(26);
			System.out.println("After filter one we have: " + lib.getAllSpectrums().size() + " candidates");
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<Spectrum> candidates = filteredLib.getSpectrumByPeaks(queryPeaks, 2);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidates, scorer1, scorer1);
			TheoreticalSpectrum t = (TheoreticalSpectrum)searcher.topSpectrum(s);
			List<Peak> annotated = searcher.getAnnotatedPeak(s, t);
			if(annotated.size() > 30){
				annotated = annotated.subList(annotated.size()-10, annotated.size());
			}
			queryPeaks = s.getTopPeaks(25, annotated);
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<Spectrum> candidates2 = filteredLib.getSpectrumByPeaks(queryPeaks, 2);
			candidates.addAll(candidates2);
			String[] peps = s.peptide.split(" & ");
			boolean passedFilter = checkPassFilter(peps[0], candidates);
			boolean passedFilter2 = checkPassFilter(peps[1], candidates);
			System.out.println("After filter two we have: " + candidates.size() + " candidates ");
			System.out.println("After filter correct peptide is retained?: " + passedFilter + "\t" + passedFilter2);			
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void testAnalyzeMatchedPeaks(){
		String spectrumFile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib.txt";
		String mixturefile = "..\\mixture_linked\\exp_sim_mixtures.id"; 
//		SpectrumLib mixlib = lib2.createMix(mixturefile, 110, 1, 0.2, 0.0000000, 1.0, 3, false);
		SpectrumLib mixlib = lib2.createRandomMix(100, 1, 0.2, 0.0000000, 1.0, 3, false);
		List<Spectrum>  specList = mixlib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		factory.setMinCharge(2);
		factory.setMaxCharge(3);
		System.out.println("Total memory used to index library: " + (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()));
		long start = (new GregorianCalendar()).getTimeInMillis();
		System.out.println("mixlib has size: " + specList.size());
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(specList, null, null);
			String[] peps = s.peptide.split(" & ");
			TheoreticalSpectrum t = new TheoreticalSpectrum(peps[0], new String[]{"b"}, new String[]{"y"});
			List<Peak> annon = searcher.getAnnotatedPeak(s, t);
			//s.getPeak().removeAll(annon);
			//s.computePeakRank();
			TheoreticalSpectrum t2 = new TheoreticalSpectrum(peps[1]);
			t2.analyzeAnnotation(s, peps[1]);
			//t.analyzeMixtureAnnotation(s, peps[0], peps[1]);
		}
	}
	
	public static void testAnalyzeMatchLinkedPeaks(){
		String spectrumFile2 = "..\\mixture_linked\\linked_peptide_spectra2.mgf";
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		List<Spectrum>  specList = lib2.getSpectrumList();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.getPeak().size() < 10){
				continue;
			}
			//s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			//SpectrumLibSearcher searcher = new SpectrumLibSearcher(specList, null, null);
			String[] peps = s.peptide.split(" & ");
			//TheoreticalSpectrum t = new TheoreticalSpectrum(peps[0], new String[]{"b"}, new String[]{"y"});
			//List<Peak> annon = searcher.getAnnotatedPeak(s, t);
			//s.getPeak().removeAll(annon);
			//s.computePeakRank();
			TheoreticalSpectrum t2 
				= TheoreticalSpectrum.getLinkedTheoreticalSpectrum(peps[0], peps[1], 
						(short)s.charge, 4, 4);
			//t2.analyzeAnnotation(s, peps[1]);
			t2.analyzeMixtureAnnotation(s, peps[0], peps[1]);
		}
	}
	
	public static void main(String[] args){
		//testLookUpSpectrum();
		//testLookUpMixtureSpectrum();
		//testIterativeLookUpMixtureSpectrum();
		//testAnalyzeMatchedPeaks();
		testAnalyzeMatchLinkedPeaks();
	}
}
