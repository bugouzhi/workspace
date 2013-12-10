package org.Spectrums;

public class SpectrumPepFilter implements SpectrumFilter{
	String peptidePattern = "";
	public SpectrumPepFilter(String pattern){
		this.peptidePattern = pattern;
	}
	
	@Override
	public boolean accept(Spectrum s) {
		//System.out.println(s.peptide);
		//System.out.println("pattern: " + this.peptidePattern);
		return s.peptide.matches(this.peptidePattern);
	}

}
