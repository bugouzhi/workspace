package org.Spectrums;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.systemsbiology.jrap.stax.Scan;

/**
 * Some utility function for SWATH data
 * @author Jian Wang
 *
 */
public class SWATHUtils {
	public static final int DALTON = 1;
	public static final int PPM = 2;
	
	
	public static boolean checkPrecursor(double ms2, double precursor, double windowWidth){
		 double diff = precursor - ms2 + 4.0;
		 return diff > 0 && diff < windowWidth;
	}
	
	public static boolean checkMass(double ms2, double precursor, double tolerance, int mode){
		return SpectrumUtil.checkMass(ms2, precursor, tolerance, mode);
	}
	
	
	public static int[] peakIntDistr(Spectrum s, double min, double max, double binWidth){
		List<Peak> sorted = new ArrayList();
		sorted.addAll(s.getPeak());
		Collections.sort(sorted, PeakIntensityComparator.comparator);
		int bins = (int)Math.ceil((max - min)/binWidth);
		int[] counts = new int[bins];
		for(int i = 0; i < sorted.size(); i++){
			
		}
		return null;
	}
	
	public static double getRT(Scan s){
		String rt = s.getHeader().getRetentionTime();
		return Double.parseDouble(rt.substring(2, rt.length()-1));
	}
	
	public static int getSWATHMS1Scan(Spectrum s){
		return (int)(s.scanNumber / 35.0)*35+1;
	}
}

