package MSPLIT;

public class MixtureSpectrum extends SimpleSpectrum{
	private int[] charges;
	private double[] parentMasses;
	private int precursorCount;
	
	//generate mixture spectrum from two simple spectrum
	public MixtureSpectrum(Spectrum spectrum1, Spectrum spectrum2, double scale1, double scale2){
		//for simplicity we set all the spectrum attribute to those of the first spectrum
		//in the mixture, presumably this correspond to the highest abundant peptides in mixture spectrum
		SimpleSpectrum s1 = (SimpleSpectrum)spectrum1;
		SimpleSpectrum s2 = (SimpleSpectrum)spectrum2;
		this.setCharge(s1.getCharge());
		this.setParentMass(s1.getParentMass());
		this.setSpectrumName(s1.getSpectrumName());
		this.charges = new int[]{s1.getCharge(), s2.getCharge()};
		this.parentMasses = new double[]{s1.getParentMass(), s2.getParentMass()};
		this.precursorCount = 2;
		this.setSpectrumName(s1.getSpectrumName() + "\t" + s2.getSpectrumName());
		Peak p1, p2;
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < s1.peaks.size() && j < s2.peaks.size()){
			p1 = s1.peaks.get(i);
			p2 = s2.peaks.get(j);
			mz1 = p1.getMass();
			mz2 = p2.getMass();
			if(mz1 < mz2){
				this.peaks.add(new SimplePeak(mz1, p1.getIntensity()*scale1));
				i++;
			}else if(mz1 == mz2){
				this.peaks.add(new SimplePeak(mz1, p1.getIntensity()*scale1 + p2.getIntensity()*scale2));
				i++;
				j++;
			}else {
				this.peaks.add(new SimplePeak(mz2, p2.getIntensity()*scale2));
				j++;
			}
		}
		//appending any remaining peaks 
		while(i < s1.peaks.size()){
			p1 = s1.peaks.get(i);
			this.peaks.add(new SimplePeak(p1.getMass(), p1.getIntensity()*scale1));
			i++;
		}
		while(j < s2.peaks.size()){
			p2 = s2.peaks.get(j);
			this.peaks.add(new SimplePeak(p2.getMass(), p2.getIntensity()*scale2));
			j++;
		}
	}

	public int[] getCharges() {
		return charges;
	}

	public void setCharges(int[] charges) {
		this.charges = charges;
	}

	public double[] getParentMasses() {
		return parentMasses;
	}

	public void setParentMass(double[] parentMass) {
		this.parentMasses = parentMass;
	}

	public int getPrecursorCount() {
		return precursorCount;
	}

	public void setPrecursorCount(int precursorCount) {
		this.precursorCount = precursorCount;
	}
}	
