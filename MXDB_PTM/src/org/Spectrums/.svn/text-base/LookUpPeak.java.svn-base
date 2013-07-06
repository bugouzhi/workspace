package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

/**
 * Compute look up peaks of a give spectrum, use to filter
 * out candidates
 * @author Jian Wang
 *
 */
public class LookUpPeak {
	 private Spectrum s;
	 public LookUpPeak(Spectrum s){
		 this.s = s;
	 }
	 
	 /**
	  * Find complementary peaks pairs that sums up to parent mass
	  * @return
	  */
	 public List<Pair> getComplementaryPeaks(double massTolerance){
		 List<Pair> peakPairs = new ArrayList();
		 List<Peak> peaks = this.s.getPeaks();
		 int i = 0;
		 int j = peaks.size()-1;
		 Peak p1, p2;
		 double pmMass = s.parentMass * s.charge - s.charge*Mass.PROTON_MASS;
		 while(i < j  && i < peaks.size() && j > 0){
			 p1 = peaks.get(i);
			 p2 = peaks.get(j);
			 //if total mass get too small, advance i and rewind j
			 if(p1.getMass() + p2.getMass() < pmMass){
				i++;
				do{
					p1 = peaks.get(i);
					p2 = peaks.get(j);
					j++;
				}while(j < peaks.size() && p1.getMass() + p2.getMass() < pmMass);
			 }
			 
			 //check if pair is  within tolerance, 
			 if(Math.abs(p1.getMass() + p2.getMass() - pmMass) < 2*massTolerance){
				 Pair p = new Pair(p1, p2);
				 peakPairs.add(p);
			 }
			 //try another pair
			 j--;
		 }
		 return peakPairs;
	 }
	 
	 /**
	  * Find peaks that are an amino acid apart by mass
	  * @param massTolerance
	  * @return
	  */
	 public List<Pair> getConsecutivePeaks(double massTolerance){
		 List<Pair> peakPairs = new ArrayList();
		 List<Peak> peaks = this.s.getPeaks();
		 for(int i = 0; i < peaks.size(); i++){
			 Peak p1 = peaks.get(i);
			 int j = i+1;
			 Peak p2 = peaks.get(j);
			 while(j < peaks.size() && p2.getMass() - p1.getMass() < Mass.maxAAMass){
				 if(Mass.isAAMass(p2.getMass() - p1.getMass())){
					 Pair p = new Pair(p1, p2);
				 }
				 j++;
			 }
			 
		 }
		 return peakPairs;
	 }
	 
	 public void computePeakScore(){
		 
	 }
	 
	 
	 
	 

}
