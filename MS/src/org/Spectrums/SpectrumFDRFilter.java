package org.Spectrums;

public class SpectrumFDRFilter implements SpectrumFilter{
	public static int PSM = 0;
	public static int PEP = 1;
	double minFDR = 0.0;
	double maxFDR = 0.1;
	int mode = PSM;
	
	public  SpectrumFDRFilter (double max, int mode){
		this.maxFDR = max;
		this.mode = mode;
	}
	
	public  SpectrumFDRFilter (double max, double min,  int mode){
		this.minFDR = min;
		this.maxFDR = max;
		this.mode = mode;
	}
	@Override
	public boolean accept(Spectrum s) {
		AnnotatedSpectrum as = (AnnotatedSpectrum)s;
		double fdr = 0;
		if(mode == PSM)
			fdr = (Double)as.getAnnotation().get("fdr");
		if(mode == PEP)
			fdr = (Double)as.getAnnotation().get("pepfdr");
		//System.out.println("fdr: " + fdr);
		return fdr <= this.maxFDR && fdr >= this.minFDR;
	}

}
