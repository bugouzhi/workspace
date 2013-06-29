package org.Spectrums;

import java.util.Set;
import java.util.Iterator;

public class MixtureSpectrumScorer extends SimpleProbabilisticScorer{
	public static boolean detail = false;
	public MixtureSpectrumScorer(PeakComparator comp) {
		super(comp);
	}
	
	@Override
	public double compare(Spectrum s1, Spectrum s2) {
		if(!(s1 instanceof TheoreticalSpectrum)){
			throw new IllegalArgumentException("First argument must be a theoretical spectrum");
		}
		TheoreticalSpectrum t = (TheoreticalSpectrum) s1;
		//System.out.println("comparing: " + t.peptide +"\t" + s2.peptide + "linked position ");
		SimpleMatchingGraph g = t.getMatchGraph(s2, this.matchTolerance);
		t.refineMatchedSpectrum(g, s2);
		double score1 = this.computeScore(g, false, this.includeNoise);
		reversePeptide(g);
//		double score2 = score1;
		double score2 = this.computeScore(g, false, this.includeNoise);
		reversePeptide(g);
		//System.out.println("score1: " + score1 + " score2: " + score2);
		if(score1 > score2){
			return score1;
		}else{
			return score2;
		}
	}
	
	public void reversePeptide(SimpleMatchingGraph g){
		Set nodes = g.vertexSet(2);
		for(Iterator it = nodes.iterator(); it.hasNext();){
			MixturePeak m = (MixturePeak)it.next();
			if(m.getPeptideIndex() == 0){
				m.setPeptideIndex(1);
			}else{
				m.setPeptideIndex(0);
			}
		}
	}
	

}
