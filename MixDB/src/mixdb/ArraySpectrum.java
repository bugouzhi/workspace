package mixdb;

import org.Spectrums.Peak;
import org.Spectrums.Spectrum;

/**
 * An array-based implementation of spectrum
 * allow fast acces to peaks
 * @author Jian Wang
 *
 */
import org.Spectrums.ArrayUtils;
/**
 * An implementation of the spectrum back by an array.  The spectrum
 * is basically a arrays of mass and intensity information
 * This allow for much faster processing and analyzing the spectrum compared to
 * the Peak-objected based spectrum
 * @author Jian
 *
 */
public class ArraySpectrum extends Spectrum{
	public static int MASS = 0;
	public static int INTENSITY = 1;
	private double[][] massIntensityList=null;
	
	
	public ArraySpectrum(){
		
	}
	
	public ArraySpectrum(double[][] massIntensityList, double parentmass, int charge){
		this.massIntensityList = massIntensityList;
		this.parentMass = parentmass;
		this.setCharge(charge);
	}
	
	public ArraySpectrum(Spectrum s){
		super(s);
	}
	
	public double[][] getMassIntensityList() {
		return massIntensityList;
	}
	public void setMassIntensityList(double[][] massIntensityList) {
		this.massIntensityList = massIntensityList;
	}
	
	public int getCharge() {
		return charge;
	}
	public void setCharge(int charge) {
		this.charge = charge;
	}
	
	public static ArraySpectrum getArraySpectrum(Spectrum s){
		double[][] massIntList = new double[2][s.getPeak().size()];
		for(int i = 0; i < s.getPeak().size(); i++){
			Peak p = s.getPeak().get(i);
			massIntList[ArraySpectrum.MASS][i] = p.getMass();
			massIntList[ArraySpectrum.INTENSITY][i] = p.getIntensity();
		}
		return new ArraySpectrum(massIntList, s.parentMass, s.charge);
	}
	
	public static ArraySpectrum getRankSpectrum(Spectrum s){
		ArraySpectrum acopy = new ArraySpectrum(s);
		int[] rankInterval = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30, 35, 
				40, 45, 50, 55,60,65, 70,75,80, 85, 90,100,110, 120, 130, 140, 150,	2000};
		double[][] massIntList = new double[2][s.getPeak().size()];
		for(int i = 0; i < s.getPeak().size(); i++){
			Peak p = s.getPeak().get(i);
			massIntList[ArraySpectrum.MASS][i] = p.getMass();
			massIntList[ArraySpectrum.INTENSITY][i] = ArrayUtils.getIntervalIndex(p.getRank(), rankInterval)+1;
		}
		acopy.setMassIntensityList(massIntList);
		return acopy;
	}
	
}
