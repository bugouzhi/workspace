package MSPLIT;

/**
 * A spectrum is a list of peaks with further requirement
 * that peaks are ordered according to their mass
 */

import java.util.Iterator;
import java.util.List;

public interface Spectrum{
	/**
	 * @return peaks according to their mass
	 */
	Iterator<Peak> peakIterator();
		                          
	/**                       
	 * @param position	
	 * @return the peak at position
	 */
	Peak getPeak(int position);
	
	/**
	 * add a peak
	 */
	void addPeak(Peak p);
	
	/**
	 * total number of peaks in the spectrum
	 */
	int numOfPeaks();
	
	/**
	 * remove a peak
	 */
	Peak removePeak(int i);
	
	boolean removePeak(Peak p);
	
	public int getCharge();
	
	public double getParentMass();
	
	public String getSpectrumName();
	
	public void setParentMass(double parentMass);
	
	public void setCharge(int charge);
	
	public void setSpectrumName(String spectrumName);
	public void setPeaks(List<Peak> peakList);
}
