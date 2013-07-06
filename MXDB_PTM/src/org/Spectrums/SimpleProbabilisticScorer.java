package org.Spectrums;

import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.Set;
import java.util.Map;
import java.util.List;
import java.util.HashMap;

import mixdb.ArraySpectrum;
import mixdb.ArrayTheoreticalSpectrum;
public class SimpleProbabilisticScorer implements SpectrumComparator{
	protected PeakComparator comp;
	protected boolean includeNoise = true;
	protected double matchTolerance = 0.5;
	protected int minMatchedPeak = 1;
	public int matchToleranceMode = 0;
	/**
	 * 
	 * @param comp
	 */
	public SimpleProbabilisticScorer(PeakComparator comp){
		this.comp = comp;
	}
	
	public double computeScore(SimpleMatchingGraph g, boolean detail, boolean includeNoise){
		double matchScore = 0.0, unMatchScore = 0.0, noiseScore = 0.0;
		int matchCount = 0, unMatchCount = 0, noiseCount = 0;
		int matchedPeak = 0;
		Map<Peak, Double> scoreTable = new HashMap<Peak, Double>();
		for(Iterator<? extends Peak> it = g.vertexSet(1).iterator(); it.hasNext();){
			Peak experimental = it.next();
			Set<? extends Peak> neighbors = g.getNeighborSet(experimental, 1);
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
			Set<? extends Peak> neighbors = g.getNeighborSet(theoretical, 2);
			if(neighbors.size() == 0){
				if(((LabelledPeak)theoretical).getPep().getPeptide().contains("Z")){
					continue;
				}
				unMatchScore += this.computeNotMatchScore(theoretical);
				unMatchCount++;
			}
		}
		if(detail){
			System.out.println("matching score breakdown: " + matchScore + "\t" + unMatchScore+"\t"+noiseScore +"\t" + (matchScore+unMatchScore));
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
		SimpleMatchingGraph g = null;
//		return fastCompare(s1, s2, false);
		int oldCharge = t.charge;
		int newCharge = oldCharge > 4 ? 4 : oldCharge;
		//System.out.println("Tolerance mode: " + this.matchToleranceMode + " tolerance: " + this.matchTolerance);
		if(s1 instanceof ArraySpectrum){
			g = t.getMatchGraph(s2, this.matchTolerance, false, this.minMatchedPeak);
		}else if(matchToleranceMode == 0){
			g = t.getMatchGraph(s2, this.matchTolerance, false, this.minMatchedPeak);
		}else {
			g = t.getPPMMatchGraph(s2, this.matchTolerance, false, this.minMatchedPeak);
		}
		//t.matchSpectrum(s2, this.matchTolerance);
		//g = t.refineMatchedSpectrum(g, s2);
		//t.charge = oldCharge;
		return this.computeScore(g, false, this.includeNoise);
		//return Math.random();
	}
	
	//since this method is called so often, reimplement this wihtout the match graph construct
	public double fastCompare(Spectrum s1, Spectrum s2, boolean DEBUG) {
		if(!(s1 instanceof TheoreticalSpectrum)){
			throw new IllegalArgumentException("spectrum 1 being compared must be of type ArrayTheoreticalSpectrum");
		}
		
		if(!(s2 instanceof Spectrum)){
			throw new IllegalArgumentException("spectrum 2 being compared must be of type ArraySpectrum");
		}
		double score = 0.0, matchesScore = 0.0, unMatchesScore = 0.0;
		TheoreticalSpectrum spec1 = (TheoreticalSpectrum)s1;
		Spectrum spec2 = s2;
		List<Peak> pList1 = spec1.getPeak();
		List<Peak> pList2 = spec2.getPeak();
		boolean[] isMatched = new boolean[pList1.size()];
		int matchedCount = 0, matchedCount2=0, unMatchedCount = 0;
		int i = 0, j = 0;
		while(i < pList1.size() && j < pList2.size()){
			//System.out.println(i + "\t" + j);
			Peak p1 = pList1.get(i);
			Peak p2 = pList2.get(j);
			double mass2 = p2.getMass();
			double mass1 = p1.getMass();
			if(mass2 - mass1 > matchTolerance){
				i++;
				continue;
			}
			if(mass1 - mass2 > matchTolerance){
				j++;
				continue;
			}
			double peakMaxScore = -10000.0;
			double mass = mass1;
//			for(int p = i; mass - mass2 < this.massTolerance && p < massIntList1[0].length; 
//				p++, mass = massIntList1[ArraySpectrum.MASS][p]){
			int p = i;
			while(p < pList1.size()){
				Peak p11 = pList1.get(p);
				mass = p11.getMass();
				if(mass - mass2 > matchTolerance){
					break;
				}
				isMatched[p]=true;
				double currScore = 	this.comp.compare(p11, p2);
				//double errorScore = getErrorScore(massIntList1[ArraySpectrum.MASS][p]-massIntList2[ArraySpectrum.MASS][j],
				//		(int)massIntList2[ArraySpectrum.INTENSITY][j]);
				//currScore += errorScore;
				peakMaxScore = peakMaxScore > currScore ? peakMaxScore : currScore;
				if(DEBUG){
					System.out.println("peak: " + pList2.get(j) + "\t" + pList1.get(p)+ "\tcurrent score: " + currScore);// + errorScore);
				}
				p++;
			}
			
			if(peakMaxScore != -10000){
				matchedCount++;
			}
			matchesScore += peakMaxScore;
			if(DEBUG){
				System.out.println("current max: " + peakMaxScore);
			}
			j++;
			continue;
		}
		//take care of unmatched theoretical fragments
		for(int k = 0; k < isMatched.length; k++){
			if(!isMatched[k]){
				LabelledPeak p = (LabelledPeak)pList1.get(k);
				if(p.getPep().getPeptide().contains("Z")){    //don't want to penalize for dummy substrate peptides
					continue;
				}
				unMatchesScore += this.comp.compare(pList1.get(k), null);  //zero means unmatched peaks
				unMatchedCount++;
			}else{
				matchedCount2++;
			}
		}
		score = matchesScore + unMatchesScore;
		if(DEBUG){
			System.out.println("score: " + score + "\t" + matchesScore + "\t" + unMatchesScore);
			System.out.println("matched count: " + matchedCount + "(" + matchedCount2 + ")" + "\t" + unMatchedCount);
		}
		return score;
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
		SpectrumIonRankLearner learner = new SpectrumIonRankLearner(lib1);
		PeakComparator peakscorer = learner.createComparatorSet();
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < 5000; i++){
			Spectrum s = specList.get(i);
			if(s.charge >=2 && s.charge < 4 || true){
				TheoreticalSpectrum t = new TheoreticalSpectrum(s.peptide);
				System.out.println("Spectrum: " + s.peptide + " has score: " + scorer.compare(t, s));
				StringBuffer buff = new StringBuffer(t.peptide.split("\\.")[0]);
				t = new TheoreticalSpectrum(buff.reverse().toString() + ".2");
				System.out.println("Spectrum: " + s.peptide + " has score: " + scorer.compare(t, s));
			}
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
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
		testSimpleScorer();
		//testRankBaseScorer();
	}
}
