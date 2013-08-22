package org.Spectrums;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * An wrapper of spectrum that enable quick query of peaks
 * in spectrum by mass or intensity
 * @author Jian
 *
 */
public class SpectrumMap{
	private Spectrum s;
	TreeMap<Double, Peak> intensityMap = new TreeMap<Double, Peak>();
	TreeMap<Double, Peak> massMap = new TreeMap<Double, Peak>();
	
	public SpectrumMap(Spectrum s){
		List<Peak> newPeakList = new ArrayList<Peak>();
		for(int i = 0; i < s.getPeaks().size(); i++){
			Peak current = s.getPeaks().get(i);
			current.setIntensity(current.getIntensity()+Math.random()*0.000000001);//in case exact same intensity
			current.setMoz(current.getMass()+Math.random()*0.000000001); //in case exact same mass
			intensityMap.put(current.getIntensity(), current); 
			massMap.put(current.getMass(), current); //in case exact same masses
		}
	}
	
	public boolean checkPeak(double mass, double tolerance){
		//System.out.println("checking MS1");
		SortedMap<Double, Peak> foundPeaks = massMap.subMap(mass - tolerance, mass + tolerance);
		//System.out.println("found: " + foundPeaks.size());
		return  foundPeaks.size() > 0;
	}
	
	public boolean checkChargedPeaks(double mass, double tolerance, int charge){
		//System.out.println("checking MS1");
		double c13 = Mass.C13 - Mass.C12;
		SortedMap<Double, Peak> foundPeaks = massMap.subMap(mass - tolerance, mass + tolerance);
		double isoPeak = mass + c13/charge;
		//System.out.println(isoPeak + "\t" + mass + "\t" + charge);
		SortedMap<Double, Peak> foundPeaks2 = massMap.subMap(isoPeak - tolerance, isoPeak + tolerance);
		//System.out.println("found: " + foundPeaks.size());
		return  foundPeaks.size() > 0 && foundPeaks2.size() > 0;
	}
	
	
//	public boolean checkPeaks(double[] masses, double tolerance){
//		SortedMap<Double, Peak> foundPeaks = massMap.subMap(masses[0] - tolerance, masses[masses.length] + tolerance);
//		//int i 
//		for(Iterator<Peak> it = foundPeaks.values().iterator(); it.hasNext();){
//			Peak p = it.next();
//			
//		}
//	}
}
