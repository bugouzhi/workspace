package org.Spectrums;

import java.util.HashMap;
import java.util.Map;

import Utils.FileIOUtils;

/**
 * A simple peak comparator that score a given theoretical peak against a experimental peaks.
 * In this simple model, only the ion type of the theoretical peak matter.
 * @author Jian Wang
 *
 */
public class SimplePeakComparator implements PeakComparator{
	private Map<String, Double> matchModel;
	private Map<String, Double> unMatchModel;
	private double noiseProb = 0.05;
	private double[] noiseModel;

	//Note: we leave the control of whehter to consider peak at the spectrum level rather than here
	/**
	 * Create a scorer from parameter files
	 * @param scoreModel - this file should contain probability for present/absent of each
	 * ion type, in format of ["ion-type" + "\t" + "probability"] format. Type should be in the format 
	 * b-H20@3@3: three fields separate by @ first field is the ion type, follow by ion charge 
	 * and the peptide charge 
	 * @param noiseModel - this file contain probability for noise peaks showing up at particular 
	 * peak rank. Format is ["rank" + "\t" + "\t" + "probability"] format
	 */
	public SimplePeakComparator(String scoreModel, String noiseModel){
		this.matchModel = new HashMap<String, Double>();
		this.unMatchModel = new HashMap<String, Double>();
		createScoreModel(scoreModel);
		createNoiseModel(noiseModel);
	}
	
	private void createScoreModel(String score){
		Map<String, String> table = FileIOUtils.createTableFromFile(score, 0, 1);
		double noiseComplement = 1-this.noiseProb;
		for(java.util.Iterator<String> keys = table.keySet().iterator(); keys.hasNext();){
			String key = keys.next();
			double value = Double.parseDouble(table.get(key));
			this.matchModel.put(key, new Double(Math.log(value/this.noiseProb)));
			this.unMatchModel.put(key, new Double(Math.log((1-value)/noiseComplement)));
		}
	}
	
	private void createNoiseModel(String noiseModel){
		Map<String, String> table = FileIOUtils.createTableFromFile(noiseModel, 0, 1);
		this.noiseModel = new double[table.keySet().size()]; 
		for(java.util.Iterator<String> keys = table.keySet().iterator(); keys.hasNext();){
			String key = keys.next();		
			double value = Double.parseDouble(table.get(key));
			this.noiseModel[Integer.parseInt(key)-1] = Math.log(value);
		}
	}
	
	//we also  handle cases where the peaks is from a linked peptide
	//in that case we choose the appropriate model base on the linked charge
	//and the charge on the original peptide
	public double compare(Peak p1, Peak p2) {
		if(p1 == null){
			//System.out.println("computing noise with :" + p2.getRank());
			if(p2.getRank() < noiseModel.length){
				return noiseModel[p2.getRank()];
			}else{
				return -0.5;
			}
		}
		LabelledPeak lp = (LabelledPeak)p1;
		
		int pepCharge = lp.getPep().getCharge();
		int peakCharge = lp.getCharge();
		if(peakCharge > pepCharge){
			peakCharge = peakCharge - pepCharge - 1;
		}
		if(peakCharge == 0){
			return 0;
		}
		if(p2 == null){
			//if(p1.getMass() < 200 || p1.getMass() > 1600){   //if peak out of range, we don't really expect it to appear
			//	return 0.0;
			//}
//			System.out.println(lp.getType());
			return this.unMatchModel.get(lp.getType() + "@" + peakCharge + 
				"@" + pepCharge).doubleValue();
		}else{
//			System.out.println(lp.getType() + "@" + peakCharge + 
//					"@" + pepCharge);
			return this.matchModel.get(lp.getType() + "@" + peakCharge + 
					"@" + pepCharge).doubleValue();
		}
	}
	
	
}
