package Spectrum;
/**
 * Array implementation of a spectrum
 * It is used for fast processing of spectrum in constrat to looping through a list of peak
 * objects in the classic spectrum object (Note: only computing cosine is supported, add in other methods as needed)
 * @author Jian
 *
 */
public class ArraySpectrum extends Spectrum{
	float[][] peaks;
	public ArraySpectrum(Spectrum s){
		this.charge = s.charge;
		this.parentMass = s.parentMass;
		this.peptide = s.peptide;
		this.protein = s.protein;
		this.scanNumber = s.scanNumber;
		this.score = s.score;
		this.spectrumName = s.spectrumName;
		peaks = new float[2][s.getPeak().size()];
		for(int i = 0; i < s.getPeak().size(); i++){
			Peak current = s.getPeak().get(i);
			peaks[0][i] = (float)current.getMass();
			peaks[1][i] = (float)current.getIntensity();
		}
	}
	
	public double cosineSim(ArraySpectrum s1){
		double product = 0;
		double magnitude = this.magnitude(); 
		magnitude *= s1.magnitude();
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		double mz1, mz2; 
		int i = 0, j = 0;
		int length1 = this.peaks[0].length, length2 = s1.peaks[0].length;
		while(i < length1 && j < length2){
			mz1 = this.peaks[0][i];
			mz2 = s1.peaks[0][j];
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += this.peaks[1][i]* s1.peaks[1][j];
				i++;
				j++;
			}else{
				j++;
			}
		}
		return product/magnitude;
	}
}
