//A peak object represents a peak in MS spectrum
package org.Spectrums;
import java.io.Serializable;
public class Peak implements Serializable{
	public static final long serialVersionUID = 1L;
	private double massByCharge; //m/z value
	private double intensity;
	private int rank = -1;
	private int sector = 0; //we divide spectrum into sector by mass range
	
	public int getSector() {
		return sector;
	}

	public void setSector(int sector) {
		this.sector = sector;
	}

	public Peak(double moz, double intensity){
		this.massByCharge = moz; 
		this.intensity = intensity;
	}
	
	//copy constructor
	public Peak(Peak p){
		this(p.getMass(), p.getIntensity());
		this.rank = p.getRank();
	}
	
	public String toString(){
		String mz = String.format("%1$.3f", this.massByCharge);
		String Intensity = String.format("%1$.5f", this.intensity);
		if(rank  > 0){
			return "" + mz + " " + Intensity + " rank: " + this.rank;
		}else{
			
			return "" + mz + " " + Intensity;

		}
	}
	
	public double getMass(){
		return massByCharge;
	}
	
	public double getIntensity(){
		return intensity;
	}
	
	public void setMoz(double moz){
		this.massByCharge = moz;
	}
	
	public void setIntensity(double i){
		this.intensity = i;
	}
	
	//we scale the intensity by some factor
	//useful in normalizing the spectrum and 
	//creating mix-spectra with different weight
	//for the single-peptide spectra that are used
	//to create the mix spectrum
	public void scaleIntensity(double factor){
		this.intensity = factor * this.intensity;
	}
	
	public void scaleMass(double factor){
		this.massByCharge = factor * this.massByCharge;
	}
	
	public void shiftMass(double deltaMC){
		this.massByCharge = this.massByCharge + deltaMC;
	}
	
	public void shiftMassPPM(double ppm){
		this.massByCharge = this.massByCharge + ppm*this.massByCharge/1000000;
	}
	
	public int getRank(){
		return this.rank;
	}
	
	public void copyRank(Peak p){
		this.rank = p.getRank();
	}
	
	public void setRank(int rank){
		this.rank = rank;
	}
	
	public int hashCode(){
		return (int)Math.floor((this.massByCharge*1000));
	}
	
	
}
