package Spectrum;

import java.util.Comparator;

/**
 * A peak comparator based on intensity
 * @author Jian
 *
 */
public class PeakIntensityComparator implements Comparator<Peak>{
		public static PeakIntensityComparator comparator = new PeakIntensityComparator();
		public int compare(Peak p0, Peak p1) {
			if(p0.getIntensity()> p1.getIntensity()){
				return 1;
			}else if(p0.getIntensity() == p1.getIntensity()){
				return 0;
			}else{
				return -1;
			}
		}
	
}
