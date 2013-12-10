package org.Spectrums;

public class SpectrumProteinFilter implements SpectrumFilter{
	String protPattern = "";
	public SpectrumProteinFilter(String pattern){
		this.protPattern = pattern;
	}
	
	@Override
	public boolean accept(Spectrum s) {
		return s.protein.matches(this.protPattern);
	}

}