package org.Spectrums;

import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class SimpleProbabilisticScorer2 implements SpectrumComparator{
	protected PeakComparator comp;
	protected boolean includeNoise = true;
	protected double matchTolerance = 0.5;
	protected int minMatchedPeak = 10;
	/**
	 * 
	 * @param comp
	 */
	public SimpleProbabilisticScorer2(PeakComparator comp){
		this.comp = comp;
	}
	
	public double computeScore(SimpleMatchingGraph g, boolean detail, boolean includeNoise){
		double matchScore = 0.0, unMatchScore = 0.0, noiseScore = 0.0;
		int matchCount = 0, unMatchCount = 0, noiseCount = 0;
		int matchedPeak = 0;
		Map<Peak, Double> scoreTable = new HashMap<Peak, Double>();
		for(Iterator<? extends Peak> it = g.vertexSet(1).iterator(); it.hasNext();){
			Peak experimental = it.next();
			Set<? extends Peak> neighbors = g.getNeighborSet(experimental);
			if(neighbors.size() == 0){
				noiseScore += this.computeNoiseScore(experimental, g.vertexSet(1).size());
				noiseCount++;
				if(detail){
					//System.out.println("noise peak rank: " + experimental.getRank());
				}
			}else{
				double match = this.computeMatchScore(experimental, neighbors, scoreTable); 
				//System.out.println("adding score: " + match);
				matchScore += match;
				matchCount++;
				matchedPeak+=neighbors.size();
				//System.out.println("matched " + experimental + "\t" + match);
			}
		}
		
		for(Iterator<? extends Peak> it = g.vertexSet(2).iterator(); it.hasNext();){
			Peak theoretical = it.next();
			Set<? extends Peak> neighbors = g.getNeighborSet(theoretical);
			if(neighbors.size() == 0){
				unMatchScore += this.computeNotMatchScore(theoretical);
				unMatchCount++;
			}
		}
		if(detail){
			System.out.println("matching score breakdown: " + matchScore + "\t" + unMatchScore+"\t"+noiseScore);
			System.out.println("peaks group breakdown: " +  matchCount + "(" + matchedPeak + ")" + "\t" + unMatchCount+"\t"+noiseCount);
		}
		if(includeNoise){
			return matchScore + unMatchScore + noiseScore;
		}else{
			return matchScore + unMatchScore;
		}
	}
	
	protected double computeNoiseScore(Peak p, int total){
		return this.comp.compare(null, p);
	}
	
	protected double computeMatchScore(Peak p, Set<? extends Peak> neighbors, Map<Peak, Double> scoreTable){
		double max = -100.0, current = 0.0;
		LabelledPeak matchedPeak = null;
		double offset = 0.0;
		for(Iterator<? extends Peak> it = neighbors.iterator(); it.hasNext();){
			LabelledPeak lp = (LabelledPeak)it.next();
			current = comp.compare(lp, p);
//			if(scoreTable.containsKey(lp)){
//				if(scoreTable.get(lp) > current){
//					continue;
//				}else{
//					offset = scoreTable.get(lp);
//				}
//			}
			max = current > max ? current : max;
			offset = current > max ? offset : 0.0;
			matchedPeak = current > max ? matchedPeak : lp;
		}
		scoreTable.put(matchedPeak, max);
		//System.out.println("Best matched peak has mass error: " + (matchedPeak.getMass() - p.getMass()));
		return max-offset;
	}
	
//	private double computeMatchScore(Peak p, Set<? extends Peak> neighbors){
//		Map<Peptide, Double> maxMap = new HashMap();
//		double max = -100.0, current = 0.0;
//		for(Iterator<? extends Peak> it = neighbors.iterator(); it.hasNext();){
//			LabelledPeak lp = (LabelledPeak)it.next();
//			if(maxMap.containsKey(lp.getPep())){
//				max = maxMap.get(lp.getPep());
//			}else{
//				max = -100.0;
//			}
//			current = this.matchModel.get(lp.getType() + "@" + lp.getCharge() + 
//						"@" + lp.getPep().getCharge()).doubleValue();
//			max = current > max ? current : max;
//			maxMap.put(lp.getPep(), new Double(max));
//		}
//		max = 0.0;
//		for(Iterator<Double> it = maxMap.values().iterator(); it.hasNext();){
//			max += it.next().doubleValue();
//		}
//		return max;
//	}
	
	protected double computeNotMatchScore(Peak p){
		return this.comp.compare(p, null);
	}

	@Override
	public double compare(Spectrum s1, Spectrum s2) {
		if(!(s1 instanceof TheoreticalSpectrum)){
			throw new IllegalArgumentException("First argument must be a theoretical spectrum");
		}
		TheoreticalSpectrum t = (TheoreticalSpectrum) s1;
		//System.out.println("Theoretical has size: " + t.getPeak().size());
		SimpleMatchGraphFactory.resetFactory();
		ArrayListFactory.resetFactory();
		List[] matchedPeaks = t.matchSpectrum(s2, this.matchTolerance);
		double[] scores = new double[matchedPeaks[0].size()];
		for(int i = 0; i < matchedPeaks[0].size(); i++){
			scores[i] = this.comp.compare((LabelledPeak)matchedPeaks[0].get(i), (Peak)matchedPeaks[1].get(i));
		}
		double localmax = scores[0];
		double totalscore = 0;
		for(int i = 1; i < matchedPeaks[0].size(); i++){
			if(matchedPeaks[1].get(i) != matchedPeaks[1].get(i-1)){
				totalscore += localmax;
				localmax=scores[i];
			}else{
				localmax = scores[i] > localmax ? localmax : scores[i];
			}
		}
		return totalscore;
	}
	
	
	public int getMinMatchedPeak() {
		return minMatchedPeak;
	}

	public void setMinMatchedPeak(int minMatchedPeak) {
		this.minMatchedPeak = minMatchedPeak;
	}

	public static void testSimpleScorer(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		lib1.removeModSpectra();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		String ids = "..\\mixture_linked\\ecoli_peptide.ids";
		lib1.computeRank();
		List<Spectrum> specList = SpectrumUtil.getSpectra(ids, lib1);
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		//List<Spectrum>  specList = lib1.getSpectrumList();
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < 1000; i++){
			Spectrum s = specList.get(i);
			if(s.charge > 3){
				continue;
			}
			s.windowFilterPeaks(10,25);
			TheoreticalSpectrum t = new TheoreticalSpectrum(s.peptide);
			System.out.println("Spectrum: " + s.peptide + " has score: " + scorer.compare(t, s));
			StringBuffer buff = new StringBuffer(t.peptide.split("\\.")[0]);
			t = new TheoreticalSpectrum(buff.reverse().toString() + ".2");
			System.out.println("Spectrum: " + s.peptide + " has score: " + scorer.compare(t, s));
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testRankBaseScorer(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		lib1.removeModSpectra();
		lib1.computeRank();
		List<Spectrum>  specList = lib1.getSpectrumList();
		SpectrumComparator comp = SpectrumUtil.getRankBaseScorer(spectrumFile);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.charge >=2 && s.charge < 4){
				TheoreticalSpectrum t = new TheoreticalSpectrum(s.peptide);
				//System.out.println("Spectrum: " + s.peptide + " has score: " + comp.compare(t, s));
				StringBuffer buff = new StringBuffer(t.peptide.split("\\.")[0]);
				t = new TheoreticalSpectrum(buff.reverse().toString() + ".2");
				//System.out.println("Spectrum: " + s.peptide + " has score: " + comp.compare(t, s));
			}
		}
		System.out.println("matching " + specList.size() + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public PeakComparator getComp() {
		return comp;
	}

	public void setComp(PeakComparator comp) {
		this.comp = comp;
	}

	public boolean isIncludeNoise() {
		return includeNoise;
	}

	public void setIncludeNoise(boolean includeNoise) {
		this.includeNoise = includeNoise;
	}

	public double getMatchTolerance() {
		return matchTolerance;
	}

	public void setMatchTolerance(double matchTolerance) {
		this.matchTolerance = matchTolerance;
	}

	public static void main(String[] args){
		//testSimpleScorer();
		testRankBaseScorer();
	}
}
