package org.Spectrums;
import java.util.Iterator;
import java.util.Map;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;

import Utils.ArrayUtils;

public class SpectrumIonRankLearner {
	private static final int MAXLENGTH = Integer.MAX_VALUE;
	private SpectrumLib lib;
	private int[] peptideLengthInterval; //we train peptide with different length separately
	private int[] chargeInterval;
	private PeakComparator[][] models; //the model will be represent as 
	
	public SpectrumIonRankLearner(SpectrumLib lib){
		this.lib = lib;
		peptideLengthInterval = new int[]{0, 12,MAXLENGTH}; //our default parameter
		chargeInterval = new int[] {1, 2, 3, 4, 5};
		this.models = 
			new PeakComparator[this.peptideLengthInterval.length-1][this.chargeInterval.length-1];
	}
	
	/**
	 * Learn a rank-base model base on annotated library of spectra. Learning
	 * can be done separately for peptide with different length and different
	 * charge. The user specifiy how do separate the peptides by specifying
	 * the peptide-length interval and charge interval
	 * @param lib
	 * @param peptideLengthInterval
	 * @param chargeInterval
	 */
	public SpectrumIonRankLearner(SpectrumLib lib, int[] peptideLengthInterval, 
			int[] chargeInterval){
		this.lib = lib;
		this.peptideLengthInterval = peptideLengthInterval;
		this.chargeInterval = chargeInterval;
		this.models = 
			new PeakComparator[this.peptideLengthInterval.length-1][this.chargeInterval.length-1];
	}
	
	public PeakComparatorSet createComparatorSet(){
		this.getIonsStatistics();
		return new PeakComparatorSet(this.models, this.peptideLengthInterval, this.chargeInterval) ;
	}
	
	/*
	 * Divide the spectrum library into the specified subcatergories and
	 * Learn the ion probability for each sub-categories
	 */
	private void getIonsStatistics(){
		List[][] subLists = 
			new List[this.peptideLengthInterval.length-1][this.chargeInterval.length-1];
		for(int i = 0; i < subLists.length; i++){
			for(int j = 0; j < subLists[i].length; j++){
				subLists[i][j] = new ArrayList<Spectrum>();
			}
		}
		List<Spectrum> specList = lib.getAllSpectrums();
		for(int i = 0, size = specList.size(); i < size; i++){
			Spectrum s = specList.get(i);
			try{
				int chargeIndex = this.getPeptideChargeIntervalIndex(s.charge);
				int lengthIndex = this.getPeptideLengthIntervalIndex(s.peptide.length()-2);
				subLists[lengthIndex][chargeIndex].add(s);			
			}catch(IllegalArgumentException ire){
				//do nothing, if length or charge outside of our target range, we just ignore them
			}
		}
		for(int i = 0; i < subLists.length; i++){
			for(int j = 0; j < subLists[i].length; j++){
				this.models[i][j]=getIonsStatistics((List<Spectrum>)subLists[i][j]);
			}
		}
		
		for(int i = 0; i < subLists.length; i++){
			for(int j = 0; j < subLists[i].length; j++){
				((PeakRankBaseComparator)this.models[i][j]).printTable();
			}
		}
		System.out.println();
	}
	
	private int getPeptideChargeIntervalIndex(int charge){
		return ArrayUtils.getIntervalIndex(charge, this.chargeInterval);
	}
	
	private int getPeptideLengthIntervalIndex(int length){
		return ArrayUtils.getIntervalIndex(length, this.peptideLengthInterval);
	}
	
	private PeakRankBaseComparator getIonsStatistics(List<Spectrum> spectrumList){
		PeakRankBaseComparator comp = new PeakRankBaseComparator(spectrumList.get(0).charge);
		double[][]counts = new double[comp.maxIndex1()][comp.maxIndex2()];
		double total = 0.0;
		for(Iterator<Spectrum> it = spectrumList.iterator(); it.hasNext();){
			Spectrum curr = it.next();
			TheoreticalSpectrum t = new TheoreticalSpectrum(curr.getPeptide());
			SimpleMatchingGraph g = t.getMatchGraph(curr, 0.5);
			//g = t.refineMatchedSpectrum(g, curr);
			getIonsCount(g, counts, comp);
			total += curr.parentMass;
		}
		
		//counts[0][0] = Math.floor(total - ArrayUtils.sum(counts));
		counts[0][0] = ArrayUtils.sum(counts[0])/(1-0.93);   //we cheat a little and use a generic noise model rather than learn it
		ArrayUtils.addPseudoCounts(counts);
		//compute conditional probability for prob(rank | ionType)
		for(int i = counts.length-1; i >= 0; i--){  
			ArrayUtils.normalizeArray(counts[i]);
		}
		
		//obtain the log-ratio prob(r | ion) / prob(r | noise)
		for(int i = counts.length-1; i >= 0; i--){   //the reason we iterate backward is to normalize noise last
			for(int j = 0; j < counts[i].length; j++){
				counts[i][j] /= counts[0][j];
				counts[i][j] = Math.log(counts[i][j]);
			}
		}
		comp.setProbabilityModel(counts);
		return comp;
	}
	
	private void getIonsCount(SimpleMatchingGraph g, double[][] counts, PeakRankBaseComparator comp){
		for(Iterator it = g.vertexSet(2).iterator(); it.hasNext();){
			LabelledPeak lp = (LabelledPeak)it.next();
			List<Peak> neighbors = g.getNeighbors(lp);
			if(neighbors.size() > 0){
				//counts[comp.getIndex1(lp)][comp.getIndex2(neighbors.get(0))]+=1;
			}else{
//				if(comp.getIndex1(lp) < 5){
//					System.out.println(lp.getPep().getPeptide());
//				}
				counts[comp.getIndex1(lp)][0]++;
			}
		}
		
		for(Iterator it = g.vertexSet(1).iterator(); it.hasNext();){
			Peak p = (Peak)it.next();
			List<Peak> neighbors = g.getNeighbors(p);
			if(neighbors.size() == 0){
				counts[0][comp.getIndex2(p)]++;
			}else{
				double massDiff=0, min=1000;
				Peak closest = null;
				for(Iterator<Peak> iter = neighbors.iterator(); iter.hasNext();){
					Peak neigh = iter.next();
					massDiff = Math.abs(neigh.getMass() - p.getMass());
					closest = massDiff < min ? neigh : closest;
					min = massDiff < min ? massDiff : min;
				}
				counts[comp.getIndex1((LabelledPeak)closest)][comp.getIndex2(p)]++;
			}
		}
	}
	
	public static void testRankBaseLearner(){
//		String file = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
		String file = ".\\MSPLib\\Lib\\ecoli.msp";
		String annotation = ".\\mixture_linked\\trps\\result.txt";
		SpectrumLib lib = new SpectrumLib(file, "MSP");
		//lib.windowFilterPeaks(6, 25);
		lib.removeModSpectra();
		lib.computeRank();
		SpectrumIonRankLearner learner = new SpectrumIonRankLearner(lib);
		PeakComparator peakscorer = learner.createComparatorSet();
	}
	
	public static void main(String[] args){
		testRankBaseLearner();
	}
	
	
	
	
}
