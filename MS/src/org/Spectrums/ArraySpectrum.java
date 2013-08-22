package org.Spectrums;

public class ArraySpectrum extends TheoreticalSpectrum{
	double[] masses;
	Peak[] intensities;
	Spectrum s;
	public ArraySpectrum(Spectrum s){
		super();
		this.peptide = s.peptide;
		this.parentMass = s.parentMass;
		this.scanNumber = s.scanNumber;
		this.spectrumName = s.spectrumName;
		this.setPeaks(s.getPeak());
		this.setPeaks(s.getPeak());
		if(s instanceof TheoreticalSpectrum){
			this.p = ((TheoreticalSpectrum)s).p;
		}
		this.s = s;
		int length = s.getPeak().size();
		this.masses = new double[length];
		this.intensities = new Peak[length];
		for(int i = 0; i < length; i++){
			Peak current = s.getPeak().get(i);
			this.masses[i] = current.getMass();
			this.intensities[i] = current;
		}
	}
}
