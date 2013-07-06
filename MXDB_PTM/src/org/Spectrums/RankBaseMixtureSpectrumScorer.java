package org.Spectrums;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * Rank base scorer for mixture spectrum, due to the rank-base
 * probability model is base on single peptide, ranks in mixture
 * spectrum might not be corrected, need to adjust the peak rank
 * base on the two peptides being considered
 * @author Jian Wang
 *
 */
public class RankBaseMixtureSpectrumScorer extends SimpleProbabilisticScorer{

	public RankBaseMixtureSpectrumScorer(PeakComparator comp){
		super(comp);
	}
	
	public double compare(Spectrum s1, Spectrum s2){
		if(!(s1 instanceof TheoreticalSpectrum)){
			throw new IllegalArgumentException("First argument must be a theoretical spectrum");
		}
		TheoreticalSpectrum t = (TheoreticalSpectrum) s1;
		SimpleMatchingGraph g = t.getMatchGraph(s2, this.matchTolerance);
		//g = t.refineMatchedSpectrum(g, s2);
		return this.computeScore(g, false, this.includeNoise, s1, s2);
	}
	
	public double computeScore(SimpleMatchingGraph g, boolean detail, boolean includeNoise, Spectrum s1, Spectrum s2){
		double matchScore = 0.0, unMatchScore = 0.0, noiseScore = 0.0;
		int matchCount = 0, unMatchCount = 0, noiseCount = 0;
		//we compute the score for each peptide separately
		String[] peptides = getMixturePeptides(s1);
		Map<Peak, Integer>[] newRanks = new Map[peptides.length];
		for(int i = 0; i < peptides.length; i++){
			newRanks[i]= recomputeRank(g, peptides[i]);
		}
		//newRanks[1] = recomputeRank(g, peptides[1]);
		
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
				double match = this.computeMatchScore(experimental, neighbors, peptides, newRanks); 
				//System.out.println("adding score: " + match);
				matchScore += match;
				matchCount++;
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
			System.out.println("peaks group breakdown: " + matchCount + "\t" + unMatchCount+"\t"+noiseCount);
		}
		if(includeNoise){
			return matchScore + unMatchScore + noiseScore;
		}else{
			return matchScore + unMatchScore;
		}
	}
	
	public String[] getMixturePeptides(Spectrum s){
		String[] peptides = s.getPeptide().split(" & ");
		String[] peptides2 = new String[peptides.length];
		for(int i = 0; i < peptides2.length; i++){
			peptides2[i] = peptides[i].split("\\.")[0];
		}
		return peptides2;
	}
	
	public Map<Peak, Integer> recomputeRank(SimpleMatchingGraph g, String peptide){
		List<Peak> reRankPeaks = new ArrayList();
		Map<Peak, Integer> newRanks = new HashMap<Peak, Integer>();
		for(Iterator<? extends Peak> it = g.vertexSet(1).iterator(); it.hasNext();){
			Peak experimental = it.next();
			Set<? extends Peak> neighbors = g.getNeighborSet(experimental);
			if(neighbors.size() > 0){
				for(Iterator<? extends Peak> iter = neighbors.iterator(); iter.hasNext();){
					LabelledPeak lp = (LabelledPeak)iter.next();
					if(lp.getPep().getPeptide().equals(peptide)){
						reRankPeaks.add(experimental);
						break;
					}
				}
			}else{
				reRankPeaks.add(experimental);
			}
		}
		Collections.sort(reRankPeaks, new PeakIntensityComparator());
		for(int i = 0, size = reRankPeaks.size(); i < size; i++){
			Peak p = reRankPeaks.get(i);
			newRanks.put(p, new Integer(size-i));
		}
		return newRanks;
	}
	
	protected double computeMatchScore(Peak p, Set<? extends Peak> neighbors, String[] peptides, Map<Peak, Integer>[] newRanks){
		double max = -100.0, current = 0.0;
		for(Iterator<? extends Peak> it = neighbors.iterator(); it.hasNext();){
			LabelledPeak lp = (LabelledPeak)it.next();
			boolean found = false;
			for(int i = 0; i < peptides.length; i++){
				//System.out.println("theo-peak belongs to: " + lp.getPep().getPeptide());
				if(lp.getPep().getPeptide().equals(peptides[i])){
					int oldRank = p.getRank();
					//System.out.println("old rank is: " + oldRank);
					//System.out.println("new rank is: " + newRanks[i].get(p).intValue());
					p.setRank(newRanks[i].get(p).intValue());
					current = comp.compare(lp, p);
					max = current > max ? current : max;
					p.setRank(oldRank);
					found = true;
				}
			}
			if(!found){
				System.err.println("warning: " + lp.getPep().getPeptide());
			}
		}
		return max;
	}
}
