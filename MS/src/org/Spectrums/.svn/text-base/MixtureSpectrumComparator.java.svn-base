package org.Spectrums;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

//compare theoretical mixture spectrum generated from library
//and observed mixture spectrum, try to understand their similarity and difference
public class MixtureSpectrumComparator {
	private Spectrum observed;
	private Spectrum simulated;
	public MixtureSpectrumComparator(Spectrum observed, Spectrum simulated ){
		this.observed = observed;
		this.simulated = simulated;
	}

	public void compareAnnotedSpectrum(){
		this.observed.computePeakRank();
		this.simulated.computePeakRank();
		String[] peps = this.simulated.peptide.split(" & ");
		TheoreticalSpectrum th = new TheoreticalSpectrum(peps[0], peps[1]);	
		SimpleMatchingGraph g = th.getMatchGraph(simulated, 0.5);
		List<Peak> peakList = new ArrayList<Peak>();
		peakList.addAll(g.vertexSet(1));
		Collections.sort(peakList, PeakIntensityComparator.comparator);
		SimpleMatchingGraph g2 = SpectrumUtil.constructMatchingGraph(this.observed, this.simulated, 0.5);
		for(int i = 0; i < peakList.size(); i++){
			Peak current = peakList.get(i);
			List neighbors = g.getNeighbors(current);
			if(neighbors.size() > 0){
				List neighbors2 = g2.getNeighbors(current);
				for(int j = 0; j < neighbors2.size(); j++){
					System.out.println(current + " ~ " + neighbors2.get(j));
				}
			}else{
				System.out.println(current + " ~ ");
			}
		}
	}
	
	public void compareSpectrumStatistics(){
		String[] peps = this.simulated.peptide.split(" & ");
		TheoreticalSpectrum th = new TheoreticalSpectrum(peps[0], peps[1]);
		double[] stats1 = th.analyzeMixtureAnnotation(this.observed, peps[0], peps[1]);
		double[] stats2 = th.analyzeMixtureAnnotation(this.simulated, peps[0], peps[1]);
		System.out.println("stat1: " + stats1[0] +  "\t" + stats1[2] +  "\t" + stats1[4]); //+ "\t" + stats1[6] + "\t" + stats1[8]);
		System.out.println("stat2: " + stats2[0] +  "\t" + stats2[2] + "\t" +  stats2[4]); //+ "\t" + stats2[6] + "\t" + stats2[8]);
	}
	
	
	public void getMatchCounts(){
		this.observed.computePeakRank();
		this.simulated.computePeakRank();
		String[] peps = this.simulated.peptide.split(" & ");
		TheoreticalSpectrum th = new TheoreticalSpectrum(peps[0], peps[1]);	
		getMatchCounts(this.observed, th);
		getMatchCounts(this.simulated, th);
		//getMatchedIntensity(this.observed, th);
		//getMatchedIntensity(this.simulated, th);
	}
	
	public void getMatchCounts(Spectrum s1, TheoreticalSpectrum t1){
		String[] peps = t1.peptide.split(" & ");
		peps[0] = peps[0].split("\\.")[0];
		peps[1] = peps[1].split("\\.")[0];
		SimpleMatchingGraph g = t1.getMatchGraph(s1, 0.3);
		Set<Peak> set = g.vertexSet(1);
		int count1=0, count2=0, shareCount=0;
		for(Iterator<Peak> it = set.iterator(); it.hasNext();){
			Peak actual = it.next();
			List<Peak> neigh = g.getNeighbors(actual);
			int found1 = 0, found2 = 0;
			for(int i = 0; i < neigh.size(); i++){
				LabelledPeak lp = (LabelledPeak)neigh.get(i);
				//System.out.println(lp.getPep().getPeptide());
				if(lp.getPep().getPeptide().equals(peps[0])){
					//System.out.println("peptide1 error: " + (lp.getMass() - actual.getMass()));
					System.out.println("match graph peptide1: " + actual);
					found1 = 1;
				}
				if(lp.getPep().getPeptide().equals(peps[1])){
					//System.out.println("peptide2 error: " + (lp.getMass() - actual.getMass()));
					System.out.println("match graph peptide2: " + actual);
					found2 = 1;
				}
			}
			count1 += found1;
			count2 += found2;
			shareCount += found1*found2;
		}
		System.out.println(t1.peptide + " matched peaks breakdown: " + count1 + "\t" + count2 + "\t" + shareCount + "\t" + peps[0].length() + "\t" + peps[1].length());
	}
	
	public static double[] getMatchedIntensity(Spectrum s1, TheoreticalSpectrum t1){
		String[] peps = t1.peptide.split(" & ");
		//System.out.println("peptide one is: " + peps[0]);
		//System.out.println("peptide two is: " + peps[1]);
		//peps[0] = peps[0].split("\\.")[0];
		//peps[1] = peps[1].split("\\.")[0];
		//peps[0] = peps[0].replaceAll("[0-9\\.\\+]", "");
		//peps[1] = peps[1].replaceAll("[0-9\\.\\+]", "");
		SimpleMatchingGraph g = t1.getMatchGraph(s1, 0.5);
		Set<Peak> set = g.vertexSet(1);
		double intensity1=0, intensity2=0, shareIntensity=0;
		for(Iterator<Peak> it = set.iterator(); it.hasNext();){
			Peak current = it.next();
			List<Peak> neigh = g.getNeighbors(current);
			double found1 = 0, found2 = 0;
			for(int i = 0; i < neigh.size(); i++){
				LabelledPeak lp = (LabelledPeak)neigh.get(i);
				//System.out.println("peptide " + lp.getPep().getPeptide());
				if(lp.getPep().toString().equals(peps[0])){
					found1 = current.getIntensity();
				}
				if(lp.getPep().toString().equals(peps[1])){
					found2 = current.getIntensity();
				}
			}
			intensity1 += found1;
			intensity2 += found2;
			shareIntensity += Math.pow(found1*found2, 0.5);
		}
		double total = s1.sumMagnitude();
		intensity1 /= total;
		intensity2 /= total;
		shareIntensity /= total;
		//System.out.println(t1.peptide + " matched peaks breakdown: " + intensity1 + "\t" + intensity2 + "\t" + shareIntensity);
		return new double[]{intensity1, intensity2, shareIntensity};
	}
	
	public void compareSingleAnnotedSpectrum(){
		this.observed.computePeakRank();
		this.simulated.computePeakRank();
		TheoreticalSpectrum th = new TheoreticalSpectrum(this.simulated.peptide);	
		SimpleMatchingGraph g = th.getMatchGraph(simulated, 0.5);
		List<Peak> peakList = new ArrayList<Peak>();
		peakList.addAll(g.vertexSet(1));
		Collections.sort(peakList, PeakIntensityComparator.comparator);
		SimpleMatchingGraph g2 = SpectrumUtil.constructMatchingGraph(this.observed, this.simulated, 0.5);
		for(int i = 0; i < peakList.size(); i++){
			Peak current = peakList.get(i);
			List neighbors = g.getNeighbors(current);
			if(neighbors.size() > 0){
				List neighbors2 = g2.getNeighbors(current);
				for(int j = 0; j < neighbors2.size(); j++){
					System.out.println(current + " ~ " + neighbors2.get(j));
				}
			}else{
				System.out.println(current + " ~ ");
			}
		}
	}
	
	public static void testCompareAnnotation(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
		String spectrumFile3 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		String mixFile = "..\\mixture_linked\\mixtures.ids";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		//SpectrumLib lib1 = new SpectrumLib(spectrumFile3, "MGF");
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(true);		
		lib1.removeModSpectra();
		lib1.computeRank();
		lib2.computeRank();
		SpectrumLib mixLib = lib1.createMix(mixFile, 1000, 0.000, 1.0, 300, false);
		System.out.println("Number of mixture created: " + mixLib.getAllSpectrums().size());
		String tripletFile = "..\\mixture_linked\\mixtureAnnotation7.txt";
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); 
			String[] tokens, tokens2;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				tokens = currentLine.split("\t");
				s1 = lib2.getSpectra(tokens[0]).get(0);
				double mag = s1.sumMagnitude();
				s1.scaleSpectrum(1/mag);
				System.out.println("Getting " + tokens[1]);
				if(!mixLib.getSpectrumLibrary().containsKey(tokens[1])){
					currentLine = bf.readLine();
					continue;
				}
				s2 = mixLib.getSpectra(tokens[1]).get(0);
				s1.windowFilterPeaks(12, 25);
				//s1.filterPeaks(100);
//				String[] peps = s2.getPeptide().split(" & ");
//				System.out.println("peptide is: " + s2.getPeptide());
//				Spectrum e1 = lib1.getSpectra(peps[0]).get(0);
//				Spectrum e2 = lib1.getSpectra(peps[1]).get(0);
//				e1.windowFilterPeaks(5, 50);
//				e2.windowFilterPeaks(5, 50);
//				TheoreticalSpectrum t1 = new TheoreticalSpectrum(peps[0]);
//				TheoreticalSpectrum t2 = new TheoreticalSpectrum(peps[1]);
//				double fract1 = t1.explainedPeaks2(t1, e1, 0.5);
//				double fract2 = t2.explainedPeaks2(t2, e2, 0.5);
				//System.out.println("single stat1: " + fract1); //+ "\t" + stats1[6] + "\t" + stats1[8]);
				//System.out.println("single stat2: " + fract2); //+ "\t" + stats2[6] + "\t" + stats2[8]);
				System.out.println("experimental has peaks: " + s1.getPeak().size());
				s2.computePeakRank();
				s2.windowFilterPeaks(9, 25);
				//s2.filterPeaks(100);
				System.out.println("simulated has peaks: " + s2.getPeak().size());
				s2.computePeakRank();
				MixtureSpectrumComparator comp = new MixtureSpectrumComparator(s1, s2);
				//comp.compareAnnotedSpectrum();
				comp.getMatchCounts();
				comp.compareSpectrumStatistics();
				currentLine = bf.readLine();
				///return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void testCompareSingleAnnotation(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
		String mixFile = "..\\mixture_linked\\mixtures.ids";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(true);		
		lib1.removeModSpectra();
		lib1.computeRank();
		lib2.computeRank();
		String tripletFile = "..\\mixture_linked\\spectrumAnnotation.txt";
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); 
			String[] tokens, tokens2;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				tokens = currentLine.split("\t");
				s1 = lib2.getSpectra(tokens[0]).get(0);
				double mag = s1.sumMagnitude();
				s1.scaleSpectrum(1/mag);
				System.out.println("Getting " + tokens[1]);
				s2 = lib1.getSpectra(tokens[1]).get(0);
				s1.windowFilterPeaks(10, 50);
				s1.computePeakRank();
				s2.windowFilterPeaks(10, 50);
				s2.computePeakRank();
				MixtureSpectrumComparator comp = new MixtureSpectrumComparator(s1, s2);
				comp.compareSingleAnnotedSpectrum();
				//comp.getMatchCounts();
				currentLine = bf.readLine();
				///return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void testScorer(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
		String mixFile = "..\\mixture_linked\\mixtures.ids";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(true);		
		lib1.removeModSpectra();
		lib1.computeRank();
		lib2.computeRank();
		String tripletFile = "..\\mixture_linked\\spectrumAnnotation.txt";
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); 
			String[] tokens, tokens2;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				tokens = currentLine.split("\t");
				s1 = lib2.getSpectra(tokens[0]).get(0);
				double mag = s1.sumMagnitude();
				s1.scaleSpectrum(1/mag);
				System.out.println("Getting " + tokens[1]);
				s2 = lib1.getSpectra(tokens[1]).get(0);
				s1.windowFilterPeaks(10, 50);
				s1.computePeakRank();
				s2.windowFilterPeaks(10, 50);
				s2.computePeakRank();
				th = new TheoreticalSpectrum(s2.peptide);
				System.out.println("Got score " + scorer1.compare(th, s1) + " and " + scorer1.compare(th, s2));
				currentLine = bf.readLine();
				///return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void testMassErrorModel(){
//		String spectrumFile = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
//		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		for(Iterator<Spectrum> it = lib1.getAllSpectrums().iterator(); it.hasNext();){
			Spectrum currentSpectrum = it.next();
			TheoreticalSpectrum t = new TheoreticalSpectrum(currentSpectrum.getPeptide());
			SimpleMatchingGraph g = t.getMatchGraph(currentSpectrum, 0.5);
			Set vertices = g.vertexSet(1);
			int count = 0;
			for(Iterator<Peak> pit = vertices.iterator(); pit.hasNext();){
				Peak current = pit.next();
				List<Peak> neighbors = g.getNeighbors(current);
				for(int i = 0; i < neighbors.size(); i++){
					System.out.println("error is : " + ((current.getMass() - neighbors.get(i).getMass())*1000000/neighbors.get(i).getMass()));
				}
				count++;
				if(count > 10000){
					return;
				}
			}
		}
	}
	
	public static void testSimulatedMixture(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
		String spectrumFile3 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		String mixFile = "..\\mixture_linked\\mixtures.ids";
		//SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		SpectrumLib lib1 = new SpectrumLib(spectrumFile3, "MGF");
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		SpectrumLib lib3 = new SpectrumLib(spectrumFile3, "MGF");
		lib3.windowFilterPeaks(10, 25);
		lib3.scaleSpectrumMass(0.9995);
		lib3.toNormVector(1.0, 0.5, 2000);
		lib3.normIntensity();
		lib3.toNormVector(1.0, 0.5, 2000);
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(true);		
		lib1.removeModSpectra();
		lib1.computeRank();
		lib2.computeRank();
		SpectrumLib mixLib = lib1.createMix(mixFile, 1000, 0.000, 1.0, 300, false);
		System.out.println("Number of mixture created: " + mixLib.getAllSpectrums().size());
		String tripletFile = "..\\mixture_linked\\mixtureAnnotation7.txt";
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); 
			String[] tokens, tokens2;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				tokens = currentLine.split("\t");
				s1 = lib2.getSpectra(tokens[0]).get(0);
				double mag = s1.sumMagnitude();
				s1.scaleSpectrum(1/mag);
				System.out.println("Getting " + tokens[1]);
				if(!mixLib.getSpectrumLibrary().containsKey(tokens[1])){
					currentLine = bf.readLine();
					continue;
				}
				s1.windowFilterPeaks(10, 25);
				s1 = s1.toNormVector(1.0, 0.5, 2000);
				s1.sqrtSpectrum();
				s1 = s1.toNormVector(1.0, 0.5, 2000);
				tokens2 = tokens[1].split(" & ");
				List<Spectrum> list1 = lib3.getSpectra(tokens2[0]);
				List<Spectrum> list2 = lib3.getSpectra(tokens2[1]);
				for(int i = 0 ; i < list1.size(); i++){
					Spectrum single1 = list1.get(i);
					for(int j = 0; j < list2.size(); j++){
						Spectrum single2 = list2.get(j);
						double alpha = s1.alpha(single1, single2);
						double cosine = s1.maxScore(single1, single2, Math.pow(alpha, 0.5));
						System.out.println(s1.peptide + " similarity is: " + cosine + "\talpha:\t" + alpha);
					}
				}
				
				currentLine = bf.readLine();
				///return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void testExperimentalSimulatedMixtre(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
		String spectrumFile3 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		String mixFile = "..\\mixture_linked\\exp_sim_mixtures.id";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile3, "MGF");
		SpectrumLib lib3 = new SpectrumLib(spectrumFile, "MSP");
		lib3.removeModSpectra();
		//SpectrumLib lib2 = lib3.Divide();
		SpectrumLib lib2 = lib1;
		//lib2.toNormVector(1.0, 0.5, 2000);
		//lib3.windowFilterPeaks(10, 25);
		lib3.scaleSpectrumMass(0.9995);
		//lib1.windowFilterPeaks(5, 25);
		lib1.scaleSpectrumMass(0.9995);
		SpectrumLib mixLib = lib1.createMix(mixFile, 110, 1, 1.0, 0.0000001, 1.0, 3, false);
		SpectrumLib.runSearch(lib3, mixLib, 0.5, "", 3);
		
	}
	
	public static void main(String[] args){
		testCompareAnnotation();
		//testCompareSingleAnnotation();
		//testScorer();
		//testMassErrorModel();
		//testSimulatedMixture();
		//testExperimentalSimulatedMixtre();
	}
}
