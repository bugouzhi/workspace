package org.Spectrums;

import java.util.List;

public class LazyEvaluatedSpectrum extends TheoreticalSpectrum{
	public static String[] prefixIons = {"b"};
	public static String[] suffixIons = {"y"};
	public LazyEvaluatedSpectrum(Peptide p){
		this.peptide = p +"."+p.getCharge();
		this.p = p;
		this.charge = p.getCharge();
		this.parentMass = p.getParentmass();
	}
	
	@Override
	public List<Peak> getPeak(){
		TheoreticalSpectrum t = new TheoreticalSpectrum(this.p, TheoreticalSpectrum.prefixIons, TheoreticalSpectrum.suffixIons);
		return t.getPeak();
	}
	
	public List<Peak> getPeaks(){
		return getPeak();
	}

}
