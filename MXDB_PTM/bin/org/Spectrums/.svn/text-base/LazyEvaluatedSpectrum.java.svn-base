package org.Spectrums;

import java.util.List;

public class LazyEvaluatedSpectrum extends TheoreticalSpectrum{
	public static String[] prefixIons = {"b"};
	public static String[] suffixIons = {"y"};
	public LazyEvaluatedSpectrum(Peptide p){
		this.peptide = p +"."+p.getCharge();
		this.p = p;
		this.charge = p.getCharge();
	}
	
	@Override
	public List<Peak> getPeak(){
		TheoreticalSpectrum t = new TheoreticalSpectrum(this.p, prefixIons, suffixIons);
		return t.getPeak();
	}
	
	public List<Peak> getPeaks(){
		return getPeak();
	}

}
