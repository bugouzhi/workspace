package org.Spectrums;

import java.util.Comparator;
/**
 * Compare peak in a spectrum by its mass
 * @author Jian
 *
 */
public class PeakMassComparator implements Comparator<Peak>{
	public static PeakMassComparator comparator = new PeakMassComparator();
	private double toleracne = 0.0;
	private int mode = Mass.DIFF_DA;
	public int compare(Peak p0, Peak p1) {
		//if(checkMass(p0, p1)){
		if(p0.getMass() == p1.getMass()){
			return 0;
		}else if(p0.getMass() > p1.getMass()){
			return 1;
		}else{
			return -1;
		}
	}
	
	public PeakMassComparator(){
		this(0, Mass.DIFF_DA);
	}
	
	public PeakMassComparator(double tolerance, int mode){
		this.toleracne = tolerance;
		this.mode = mode;
	}
	
	public double getLowMass(double mass){
		if(this.mode == Mass.DIFF_DA){
			return mass - this.toleracne;
		}else{
			return mass - (mass*this.toleracne/1000000);
		}
	}
	
	public double getHighMass(double mass){
		if(this.mode == Mass.DIFF_DA){
			return mass + this.toleracne;
		}else{
			return mass + (mass*this.toleracne/1000000);
		}
	}
	
	public boolean checkMass(Peak p0, Peak p1){
		return Mass.checkMass(p0.getMass(), p1.getMass(), this.toleracne, this.mode);
	}
	

	
	
	
}
