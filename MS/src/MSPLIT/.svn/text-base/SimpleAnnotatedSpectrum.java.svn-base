package MSPLIT;

import java.util.ArrayList;

/**
 * Spectrum with simple annotation
 * e.g identified petpides, proteins etc...s
 * @author bugouzhi
 *
 */
public class SimpleAnnotatedSpectrum extends SimpleSpectrum implements AnnotatedSpectrum{
	private String peptide;
	private Peptide p;
	
	public String getPeptide() {
		return peptide;
	}
	public void setPeptide(String peptide) {
		this.peptide = peptide;
	}
	public Peptide getP() {
		return p;
	}
	public void setP(Peptide p) {
		this.p = p;
	}
	
	public Object clone(){
		SimpleAnnotatedSpectrum s = new SimpleAnnotatedSpectrum();
		s.setCharge(this.getCharge());
		s.setSpectrumName(this.getSpectrumName());
		s.setParentMass(this.getParentMass());
		s.setP(this.getP());
		s.setPeptide(this.getPeptide());
		return s;
	}
	
	@Override
	public String getAnnotations() {
		return this.peptide;
	}
	
}
