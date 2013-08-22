package MSPLIT;
/**
 * spectrum back by an array to allow
 * fast computation of cosine
 * @author Jian Wang
 *
 */
public class ArraySpectrum extends SimpleSpectrum{
	double[] masses;
	double[] intensities;
	public ArraySpectrum(Spectrum s){
		super(s);
		this.masses = new double[s.numOfPeaks()];
		this.intensities = new double[s.numOfPeaks()];
		for(int i = 0; i < s.numOfPeaks(); i++){
			Peak current = s.getPeak(i);
			this.masses[i] = current.getMass();
			this.intensities[i] = current.getIntensity();
		}
	}
	
}
