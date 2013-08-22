package MSPLIT;

import java.io.Serializable;

public class SimplePeak implements Peak, Serializable{
	public static final long serialVersionUID = 1L;
	private double mass; //m/z value
	private double intensity;
	
	public SimplePeak(double mass, double intensity){
		if(mass < 0){ 
			throw new IllegalArgumentException("mass cannot be negative");
		}
		if(intensity < 0){
			throw new IllegalArgumentException("intensity cannot be negative");
		}
		this.mass = mass;
		this.intensity = intensity;
	}

	@Override
	public double getIntensity() {
		return (double)this.intensity;
	}

	@Override
	public double getMass() {
		return (double)this.mass;
	}

	@Override
	public void setIntensity(double intensity) {
		if(intensity < 0){
			throw new IllegalArgumentException("intensity cannot be negative");

		}
		this.intensity = (float)intensity;
	}

	@Override
	public void setMass(double mass) {
		if(mass < 0){
			throw new IllegalArgumentException("mass cannot be negative");
		}
		this.mass = (float)mass;
	}
	
	public String toString(){
		return "(" + this.mass + ", " + this.intensity + ")";
	}
	public static void main(String[] args){
		Peak p = new SimplePeak(1, 1);
		System.out.println("peak is: " + p);
	}
	

}
