package MSPLIT;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.List;

import MSPLIT.SpectrumIO.MSPParser;


public class CosineSpectrumComparator implements SpectrumComparator{

	public static SpectrumComparator CosineComparator = new CosineSpectrumComparator();
	public static SpectrumComparator ProjectedCosineComparator = 
			new CosineSpectrumComparator.ProjectedCosineComparator();
	private double tolerance = 0.0;
	private double parentMassTolerance = 3000.0;
	
	public double getTolerance() {
		return tolerance;
	}

	public void setTolerance(double tolerance) {
		this.tolerance = tolerance;
	}

	public double getParentMassTolerance() {
		return parentMassTolerance;
	}

	public void setParentMassTolerance(double parentMassTolerance) {
		this.parentMassTolerance = parentMassTolerance;
	}
	//non-binned version of comparator
	public static double compare(Spectrum s1, Spectrum s2, double tolerance){
		double product = 0;
		double magnitude = SpectrumUtil.magnitude(s1);
		//System.out.println("magnitude is: " + magnitude);
		magnitude *= SpectrumUtil.magnitude(s2);
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < s1.numOfPeaks() && j < s2.numOfPeaks()){
			mz1 = s1.getPeak(i).getMass();
			mz2 = s2.getPeak(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(Math.abs(mz1 - mz2) < tolerance){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += s1.getPeak(i).getIntensity()
					* s2.getPeak(j).getIntensity();
				i++;
				j++;
			}else{
				j++;
			}
		}
		return product;///magnitude;
	}
	
	public static double compare(ArraySpectrum s1, ArraySpectrum s2, double tolerance){
		double product = 0;
		//double magnitude = SpectrumUtil.magnitude(s1);
		//System.out.println("magnitude is: " + magnitude);
		//magnitude *= SpectrumUtil.magnitude(s2);
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		double mz1, mz2; 
		int i = 0, j = 0;
		int length1 = s1.masses.length;
		int length2 = s2.masses.length;
		while(i < length1 && j < length2){
			if(s1.masses[i] < s2.masses[j]){
				i++;
			}else if(Math.abs(s1.masses[i] - s2.masses[j]) < tolerance){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += s1.intensities[i]*s2.intensities[j];
				i++;
				j++;
			}else{
				j++;
			}
		}
		return product;///magnitude;
	}

	public double compare(Spectrum s1, Spectrum s2){
		if(s1 instanceof ArraySpectrum){
			return CosineSpectrumComparator.compare((ArraySpectrum)s1, (ArraySpectrum)s2, 0.5);
		}
		if(Math.abs(s1.getParentMass() - s2.getParentMass()) > this.parentMassTolerance){
			//System.out.println("not passing filters");
			return 0.0;
		}
		double product = 0;
		double magnitude = SpectrumUtil.magnitude(s1); 
		//System.out.println("magnitude is: " + magnitude);
		magnitude *= SpectrumUtil.magnitude(s2);
		//System.out.println("magnitude is: " + magnitude);
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < s1.numOfPeaks() && j < s2.numOfPeaks()){
			mz1 = s1.getPeak(i).getMass();
			mz2 = s2.getPeak(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += s1.getPeak(i).getIntensity()
					* s2.getPeak(j).getIntensity();
				i++;
				j++;
			}else{
				j++;
			}
		}
		return product/magnitude;
	}
	
	public static double alpha(Spectrum mix, Spectrum a, Spectrum b) {		
		double A, B, C, D, E ;		
		C = dot(mix,a) ;
		D = dot(a,b) ;
		E = dot(mix,b) ;
		A = dot(a,a) ;
		B = dot(b,b) ;
		double alpha = ((B*C)-(D*E) )/ ((A*E)-(C*D));
		if(alpha < 0.1 || alpha > 10){
			//return -1*alpha;
			return 0.1;
		}else{
			if(alpha > 1){
				return 1.0/alpha;
			}else{
				return alpha;
			}
		}
	}
	
	public static double dot(Spectrum s1, Spectrum s2) {
		double product = 0;
		double mz1, mz2; 
		int i = 0, j = 0;
		
		while(i < s1.numOfPeaks() && j < s2.numOfPeaks()){
			mz1 = s1.getPeak(i).getMass();
			mz2 = s1.getPeak(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("Intensity are: " + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += Math.pow(s1.getPeak(i).getIntensity(), 2)
					* Math.pow(s2.getPeak(j).getIntensity(), 2);
				i++;
				j++;
			}else{
				j++;
			}
		}
				
		return product ;
	}
	
	public static void testComparator(){
		String libFile = "..\\MSPLIB\\Lib\\human.msp";
		MSPParser parser = new MSPParser(libFile);
		List<Spectrum> specList = new ArrayList();
		for(int i = 0; i < 5000; i ++){
			specList.add(parser.next());
		}
		int counter = 0;
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s1 = specList.get(i);
			for(int j = i+1; j < specList.size(); j++){
				Spectrum s2 = specList.get(j);
				double score = CosineComparator.compare(s1, s2);
				//System.out.println(score);
				counter++;
			}
		}
		System.out.println("matching " + counter + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testArraySpectrum(){
		String libFile = "..\\MSPLIB\\Lib\\human.msp";
		MSPParser parser = new MSPParser(libFile);
		List<ArraySpectrum> specList = new ArrayList();
		for(int i = 0; i < 5000; i ++){
			specList.add(new ArraySpectrum(parser.next()));
		}
		int counter = 0;
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			ArraySpectrum s1 = specList.get(i);
			for(int j = i+1; j < specList.size(); j++){
				ArraySpectrum s2 = specList.get(j);
				double score = CosineComparator.compare(s1, s2);
				//System.out.println(score);
				counter++;
			}
		}
		System.out.println("matching " + counter + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void main(String[] args){
		//testComparator();
		testArraySpectrum();
	}
	
	static class ProjectedCosineComparator implements SpectrumComparator{
		@Override
		public double compare(Spectrum s1, Spectrum s2) {
			double product = 0;
			double magnitude = SpectrumUtil.magnitude(s1);
			//we do the dotproduct as above, but just calculates
			//the norm for s1 vector differently, we only consider
			//those values that are non-zero in this vector
			double projectedNorm = 0.0001; //very small number avoid div-by-zero error  
			double mz1, mz2; 
			int i = 0, j = 0;
			//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
			while(i < s1.numOfPeaks() && j < s2.numOfPeaks()){
				mz1 = s1.getPeak(i).getMass();
				mz2 = s2.getPeak(j).getMass();
				if(mz1 < mz2){
					i++;
				}else if(mz1 == mz2){
					//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
					product += s1.getPeak(i).getIntensity()
						* s2.getPeak(j).getIntensity();
					projectedNorm += s2.getPeak(j).getIntensity() * s2.getPeak(j).getIntensity();
					i++;
					j++;
				}else{
					j++;
				}
			}
			//System.out.println("this norm: " + this.magnitude());
			//System.out.println("mixture projected norm: " + Math.pow(projectedNorm, 0.5));
			magnitude = magnitude * Math.pow(projectedNorm, 0.5);
			return product/magnitude;
		}
	}
	
}
