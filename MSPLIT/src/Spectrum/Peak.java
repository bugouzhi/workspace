package Spectrum;

/**
 * A peak object in a mass spectrum
 */
import java.io.Serializable;
public class Peak implements Serializable{
	public static final long serialVersionUID = 1L;
	private float massByCharge; //m/z value
	private float intensity;
	//private int rank = -1;
	
	
	public Peak(double moz, double intensity){
		this.massByCharge = (float)moz; 
		this.intensity = (float)intensity;
	}
	
	//copy constructor
	public Peak(Peak p){
		this(p.getMass(), p.getIntensity());
	}
	
	public String toString(){
//		if(rank  > 0){
//			return "" + this.massByCharge + " " + this.intensity + " rank: " + this.rank;
//		}else{
			return "" + this.massByCharge + " " + this.intensity;

//		}
	}
	
	public double getMass(){
		return massByCharge;
	}
	
	public double getIntensity(){
		return intensity;
	}
	
	public void setMoz(double moz){
		this.massByCharge = (float)moz;
	}
	
	public void setIntensity(double i){
		this.intensity = (float)i;
	}
	
	//we scale the intensity by some factor
	//useful in normalizing the spectrum and 
	//creating mix-spectra with different weight
	//for the single-peptide spectra that are used
	//to create the mix spectrum
	public void scaleIntensity(double factor){
		this.intensity = (float)factor * this.intensity;
	}
	
	public void scaleMass(double factor){
		this.massByCharge = (float)factor * this.massByCharge;
	}
	
	public void shiftMass(double deltaMC){
		this.massByCharge = (float)(this.massByCharge + deltaMC);
	}
	
	public void shiftMassPPM(double ppm){
		this.massByCharge = (float)(this.massByCharge + ppm*this.massByCharge/1000000);
	}
	
	public int getRank(){
		//return this.rank;
		return (int)intensity;
	}
	
	public void copyRank(Peak p){
		//this.rank = p.getRank();
		this.setRank(p.getRank());
	}
	
	public void setRank(int rank){
		//this.rank = rank;
		this.intensity = (float)rank;
	}
	
}
