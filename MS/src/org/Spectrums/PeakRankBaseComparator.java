package org.Spectrums;

import java.util.HashMap;
import java.util.Map;

import Utils.ArrayUtils;

public class PeakRankBaseComparator implements PeakComparator{
		private double[][] scoreModel;
		private static final int MAXRANK = Integer.MAX_VALUE; 
		private int[] peakRankInterval; //peptide ranks are separate into several intervals
		private double[] massInterval; //maybe a latter add-on featuer, to-be-implemented
		private String[] ionsType;
		private int maxCharge;
		private Map<String, Integer> peakTypeMap;
		
		public PeakRankBaseComparator(){
			this(new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,60,100,MAXRANK},
					4,
					Mass.standardIonsType);
		}
		
		public PeakRankBaseComparator(int charge){
			this(new int[]{1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,60,100,200,MAXRANK},
					charge,
					Mass.standardIonsType);
		}
		
		public PeakRankBaseComparator(int[] peakRankInterval, int maxCharge, String[] ionsType){
			//should check peakRankInterval for increasing orders
			this.peakRankInterval = peakRankInterval;
			this.maxCharge = maxCharge;
			this.ionsType = ionsType;
			peakTypeMap = new HashMap();
			createIonTypeIndex();
		}
		
		private void createIonTypeIndex(){
			int index = 1;
			for(int charge = 1; charge <= maxCharge ; charge++){
				for(int type = 0; type < ionsType.length; type++){
					this.peakTypeMap.put(ionsType[type]+charge, new Integer(index++));
				}
			}
		}
		
		//we should note that 0 index is reserved for noise
		//theoretical peaks
		public int getIndex1(LabelledPeak p){
			return peakTypeMap.get(p.getType() + p.getCharge()).intValue();
		}
		
		//experimental peaks
		public int getIndex2(Peak p){
			return getPeakRankIntervalIndex(p.getRank())+1;
		}
		
		
		public int maxIndex1(){
			return peakTypeMap.size()+1;
		}
		
		public int maxIndex2(){
			return peakRankInterval.length;
		}
		
		public double compare(Peak p1, Peak p2){
			if(p1 == null){
				return 0.0;
			}
			if(p2 == null){
				return scoreModel[getIndex1((LabelledPeak)p1)][0];
			}
			//System.out.println("rank is: " + p2.getRank());
			return scoreModel[getIndex1((LabelledPeak)p1)][getIndex2(p2)];
		}
		
		private int getPeakRankIntervalIndex(int rank){
			return ArrayUtils.getIntervalIndex(rank, this.peakRankInterval);
		}
		
		public void setProbabilityModel(double[][] model){
			this.scoreModel = model;
		}
		
		public void printTable(){ 
			System.out.println("dimesnion: " + scoreModel.length + "," + scoreModel[0].length);
			ArrayUtils.printTable(this.scoreModel);
		}
}

	

