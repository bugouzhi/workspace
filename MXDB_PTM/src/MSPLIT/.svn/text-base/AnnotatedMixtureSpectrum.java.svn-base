package MSPLIT;

public class AnnotatedMixtureSpectrum extends MixtureSpectrum{
	String[] peptides;
	public AnnotatedMixtureSpectrum(AnnotatedSpectrum spectrum1, AnnotatedSpectrum spectrum2, double scale1, double scale2){
		super(spectrum1, spectrum2, scale1, scale2);
		peptides = new String[2];
		peptides[0] = spectrum1.getAnnotations();
		peptides[1] = spectrum2.getAnnotations();
	}
	public String[] getPeptides() {
		return peptides;
	}

	public void setPeptides(String[] peptides) {
		this.peptides = peptides;
	}
	
}
