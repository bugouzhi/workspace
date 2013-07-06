package MSPLIT;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

/**
 * Various utility method for spectrum 
 * @author jian wang
 *
 */
public class SpectrumUtil{
	private static double minMass = 0.5;
	private static double maxMass = 2000;
	private static double binWidth = 0.5;
	public static Spectrum toVector(Spectrum s){
		int bins = (int)((maxMass-minMass)/binWidth) + 1;
		List <Peak>  newPeaks = new ArrayList<Peak>();
		double rightBoundary = minMass + binWidth;
		double currentValue;
		double moz;
		int j = 0;
		int i = 0;
		//System.out.println("bins: " + bins);
  		for(i = 0; i < bins; i++){
			rightBoundary = minMass + i*binWidth;
			currentValue = 0;
			while(j < s.numOfPeaks()){
				moz = s.getPeak(j).getMass();
				//System.out.println("moz " + moz);
				//System.out.println("r-edge: " + rightBoundary);
				if(moz > rightBoundary){
					break;
				}else{
					currentValue += s.getPeak(j).getIntensity(); 
					j++;
				}	
			}
			if(currentValue > 0){
				//System.out.println("creating new peaks");
				newPeaks.add(new SimplePeak(i, currentValue));
			}
		}
  		s.setPeaks(newPeaks);
		//return new SimpleSpectrum(s.getSpectrumName(), s.getParentMass(), s.getCharge(), newPeaks);
  		return s;
	}
	
	public static double magnitude(Spectrum s){
		double total = 0;
		double intensity = 0;
		for(int i = 0; i < s.numOfPeaks(); i++){
			intensity = s.getPeak(i).getIntensity();
			total += intensity*intensity;
		}
		total = Math.pow(total, 0.5);
		return total;
	}
	
	public static double sumOfPeaks(Spectrum s){
		return 0.0;
	}
	
	public static void normalize(Spectrum s){
		
	}
	
	public static void normalizeByTotalIntensity(Spectrum s){
		double mag = magnitude(s);
		scaleSpectrumIntensity(s, 1/mag);
	}
	
	public static void scaleSpectrumIntensity(Spectrum s, double scale){
		for(int i = 0; i < s.numOfPeaks(); i++){
			Peak p = s.getPeak(i);
			p.setIntensity(p.getIntensity()*scale);
		}
	}
	
	public static void scaleSpectrumMass(Spectrum s, double scale){
		for(int i = 0; i < s.numOfPeaks(); i++){
			Peak p = s.getPeak(i);
			p.setMass(p.getMass()*scale);
		}
	}
	
	public static void shiftSpectrumMass(Spectrum s, double shift){
		for(int i = 0; i < s.numOfPeaks(); i++){
			Peak p = s.getPeak(i);
			p.setMass(p.getMass()+shift);
		}
	}
	
	public static void shiftSpectrumMassByPPM(Spectrum s, double ppm){
		for(int i = 0; i < s.numOfPeaks(); i++){
			Peak p = s.getPeak(i);
			double shift = p.getMass()*ppm/1000000;
			p.setMass(p.getMass()+shift);
		}
	}
		
	public static void toRelIntensity(){
		
	}
	
	public static void sqrtSpectrum(Spectrum s){
		for(int i = 0; i < s.numOfPeaks(); i++){
			Peak p = s.getPeak(i);
			p.setIntensity(Math.pow(p.getIntensity(), 0.5));
		}
	}
	public static double residual(Spectrum s1, Spectrum s2){
		Iterator<Peak> p1 = s1.peakIterator();
		Iterator<Peak> p2 = s2.peakIterator();
		Peak peak1 = null, peak2 = null;
		int exp = 4;
		double shareIntensity = 0;
		double residIntensity = 0.001; //avoid div-by-zero 
		if(p1.hasNext() && p2.hasNext()){
			peak1 = p1.next();
			peak2 = p2.next();
		}else{
			return shareIntensity/residIntensity;
		}
		double remain;
		//System.out.println("magnitude is: " + this.magnitude());
		while(p1.hasNext() && p2.hasNext()){
			if(peak1.getMass() < peak2.getMass()){
				residIntensity += Math.pow(peak1.getIntensity(),exp);
				//System.out.println("res intensity: " + peak1.getMass() + "\t" + peak1.getIntensity());
				peak1 = p1.next();
			}else if(peak1.getMass() == peak2.getMass()){
				remain = peak1.getIntensity() - peak2.getIntensity();
				if(remain < 0){
					remain = 0;
				}
				shareIntensity += Math.pow(peak1.getIntensity(),exp);
				//System.out.println("share: " + peak2.getMass() + "\t" + peak2.getIntensity());
				//System.out.println("remaining: " + peak1.getMass() + "\t" + remain);
				//residIntensity += Math.pow(remain,exp);
				peak1 = p1.next();
				peak2 = p2.next();
			}else{
				//shareIntensity += Math.pow(peak1.getIntensity(), 2);
				peak2 = p2.next();
			}
		}
		
		if(peak1.getMass() == peak2.getMass()){ //take care of last element
				shareIntensity += Math.pow(peak1.getIntensity(), exp);
				remain = peak1.getIntensity() - peak2.getIntensity();
				if(remain < 0){
					remain = 0;
				}
				//residIntensity += Math.pow(remain, 2);
		}
		if(!p1.hasNext()){
			residIntensity += Math.pow(peak1.getIntensity(), exp);
		}
		
//		if(!p2.hasNext()){
//			shareIntensity += Math.pow(peak2.getIntensity(), 2);
//			//System.out.println("share: " + peak2.getMass() + "\t" + peak2.getIntensity());
//		}
		
		while(p1.hasNext()){
			residIntensity += Math.pow(p1.next().getIntensity(), exp);
		}
		
//		while(p2.hasNext()){
//			shareIntensity += Math.pow(p2.next().getIntensity(), 2);
//		}
			
		//System.out.println("shareIntensity: " + Math.pow(shareIntensity,0.5));
//		System.out.println("s has intensity: " + s.magnitude());
		//System.out.println("resIntensity: " + Math.pow(residIntensity,0.5));
		//double alpha = Math.pow(shareIntensity/residIntensity, 0.5); //initial estimate
		double alpha = Math.pow(residIntensity, 0.5);
		//System.out.println("alpha: " + alpha);
		alpha = alpha*alpha / (1-alpha*alpha);  //reestimate by taking into account that mixture is normalized to one
		alpha =  Math.pow(alpha, 0.5); //thus each component is down-weighted slightly
		if(alpha > 1.0){
			return 1/alpha;
		}else{
			return alpha;
		}
		//return alpha;
	}
	
	public static int explainedIntensity(Spectrum s, double percent){
//		Vector<Peak> sortedPeakList = new Vector<Peak>();
//		sortedPeakList.addAll(s1.peaks);
//		Collections.sort(sortedPeakList, new peakComparator());
//		double mag = this.magnitude();
//		mag = mag*mag;
//		int i = 0;
//		Peak p;
//		double current = 0.0;
//		for(i = sortedPeakList.size()-1; i > 0;  i--){
//			current += sortedPeakList.get(i).getIntensity()
//				*sortedPeakList.get(i).getIntensity();
//			if(current / mag > percent){
//				return sortedPeakList.size() - i;
//			}
//		}
//		return sortedPeakList.size();
		return 10;
	}
	
}
