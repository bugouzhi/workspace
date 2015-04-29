package mixgf;

import java.util.HashSet;
import java.util.Set;

import org.Spectrums.Mass;
import org.Spectrums.PeakComparator;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumComparator;
import org.Spectrums.TheoreticalSpectrum;

public class PRMSpectrumComparator implements SpectrumComparator{
	protected PeakComparator comp;
	protected boolean includeNoise = true;
	protected double matchTolerance = 0.5;
	protected int minMatchedPeak = 10;
	/**
	 * 
	 * @param comp
	 */
	public PRMSpectrumComparator(){
		this.comp = null;
	}
	@Override
	public double compare(Spectrum s1, Spectrum s2) {
		if(!(s1 instanceof TheoreticalSpectrum)){
			throw new IllegalArgumentException("First argument must be a theoretical spectrum");
		}
		if(!(s2 instanceof PRMSpectrum)){
			throw new IllegalArgumentException("Second argument must be a PRM spectrum");
		}
		TheoreticalSpectrum t = (TheoreticalSpectrum) s1;
		PRMSpectrum prmSpect = (PRMSpectrum)s2;
		double[][] base = t.computeBaseMass(t.p.getPeptide(), t.p.getPos(), t.p.getPtmmasses());
		base[0][base[0].length-1]=t.parentMass*t.charge-Mass.PROTON_MASS*t.charge-Mass.WATER;
		double[] scores = prmSpect.getScoredSpectrum(t.parentMass*t.charge);
		double totalScore = 0;
		for(int i = 0; i < base[0].length; i++){
			int index = (int)Math.round((0.9995*base[0][i]));
			if(index < scores.length){
				totalScore += scores[index];
			}
		}
		return totalScore;
	}
	
	//pairs
	public double[] compare(Spectrum[] s1, Spectrum[] s2){
		TheoreticalSpectrum t1 = (TheoreticalSpectrum) s1[0];
		TheoreticalSpectrum t2 = (TheoreticalSpectrum) s1[1];
		PRMSpectrum prmSpect1 = (PRMSpectrum)s2[0];
		PRMSpectrum prmSpect2 = (PRMSpectrum)s2[1];
		double[][] base1 = t1.computeBaseMass(t1.p.getPeptide(), t1.p.getPos(), t1.p.getPtmmasses());
		double[][] base2 = t2.computeBaseMass(t2.p.getPeptide(), t2.p.getPos(), t2.p.getPtmmasses());
		base1[0][base1[0].length-1]=t1.parentMass*t1.charge-Mass.PROTON_MASS*t1.charge-Mass.WATER;
		base2[0][base2[0].length-1]=t2.parentMass*t2.charge-Mass.PROTON_MASS*t2.charge-Mass.WATER;
		double[] scores1 = prmSpect1.getScoredSpectrum(t1.parentMass*t1.charge);
		double[] scores2 = prmSpect2.getScoredSpectrum(t1.parentMass*t1.charge);
		double totalScore = 0, score1 = 0, score2 = 0;
		Set<Integer> masses1 = new HashSet<Integer>();
		for(int i = 0; i < base1[0].length; i++){
			int index = (int)Math.round(prmSpect1.getMassIndex(base1[0][i]));
			if(index < scores1.length){
				masses1.add(index);
				totalScore += scores1[index];
				score1 += scores1[index];
				//System.out.println("score: " + index + "\t" + scores1[index]);
			}
		}
		
		for(int i = 0; i < base2[0].length; i++){
			int index = (int)Math.round((prmSpect2.getMassIndex(base2[0][i])));
			if(index < scores2.length){
				if(masses1.contains(index)){
					if(scores1[index] < scores2[index]){
						totalScore += scores2[index]-scores1[index];
						//score2 += scores2[index];
					}
					//score2 += scores2[index];
				}else{
					totalScore += scores2[index];
					score2 += scores2[index];
				}
				//System.out.println("score: " + index + "\t" + scores2[index]);
			}
		}
		return new double[]{totalScore, score1, score2};
	}
	
}
