package org.Spectrums;
/**
 * A general comparator of two peaks, note that in order to
 * account for noise model in the scoring, this interface should implements compare
 * method that accept null as argument and return an appropriate score accordingly
 * @author jian wang
 *
 */


public interface PeakComparator {
	public double compare(Peak p1, Peak p2);
}
