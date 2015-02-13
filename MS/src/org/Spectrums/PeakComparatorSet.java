package org.Spectrums;

import Utils.ArrayUtils;

/**
 * A wrapper that contain multiple PeakComparator, it
 * can dispatch the appropriate comparator base-on the peptide
 * length and charge the particular peak belong to
 * @author Jian Wang
 *
 */
public class PeakComparatorSet implements PeakComparator{
	private int[] peptideLengthInterval; //we train peptide with different length separately
	private int[] chargeInterval;
	private PeakComparator[][] models; //the model will be represent as 
	public PeakComparatorSet(PeakComparator[][] comparators, int[] peptideLengthInterval, int[] chargeInterval){
		this.models = comparators;
		this.peptideLengthInterval = peptideLengthInterval;
		this.chargeInterval = chargeInterval;
	}
	
	private int getPeptideChargeIntervalIndex(int charge){
		return ArrayUtils.getIntervalIndex(charge, this.chargeInterval);
	}
	
	private int getPeptideLengthIntervalIndex(int length){
		return ArrayUtils.getIntervalIndex(length, this.peptideLengthInterval);
	}
	
	@Override
	public double compare(Peak p1, Peak p2) {
		if(p1 == null){
			return 0.0;
		}
		LabelledPeak lp = (LabelledPeak)p1;
		return models[this.getPeptideLengthIntervalIndex(lp.getPep().getPeptide().length()-2)]
		              [this.getPeptideChargeIntervalIndex(lp.getPep().getCharge())].compare(p1, p2);
	}
	
	
}
