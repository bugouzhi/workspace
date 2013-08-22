package org.Spectrums;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Compute a.a. frequency from a set of proteins
 * @author Jian Wang
 *
 */
public class AAFreq {
	String proteinFile;
	List<String> proteins;
	Map<Character, Double> aaFreq;
	public AAFreq(String filename){
		this.proteinFile = filename;
		this.proteins = Utils.FileIOUtils.createProteinsFromFasta(filename);
		computeAAFreq();
	}
	
	public Map<Character, Double> computeAAFreq(){
		this.aaFreq = new HashMap<Character, Double>();
		for(int i = 0; i < this.proteins.size(); i++){
			String protein = proteins.get(i);
			for(int j = 0; j < protein.length(); j++){
				char aa = protein.charAt(j);
				if(aaFreq.containsKey(aa)){
					double count = aaFreq.get(aa);
					count++;
					aaFreq.put(aa, count);
				}else{
					aaFreq.put(aa, 1.0);
				}
			}
		}
		long total = 0;
		for(Iterator<Character> it = this.aaFreq.keySet().iterator(); it.hasNext();){
			char aa = it.next();
			total += aaFreq.get(aa);
		}
		
		for(Iterator<Character> it = this.aaFreq.keySet().iterator(); it.hasNext();){
			char aa = it.next();
			double freq = aaFreq.get(aa)/total;
			aaFreq.put(aa, freq);
		}
		return aaFreq;
	}
	
	public Map<Double, Double> getMassFreq(){
		Map<Double, Double> massFreq = new HashMap<Double, Double>();
		for(Iterator<Character> it = this.aaFreq.keySet().iterator(); it.hasNext();){
			char aa = it.next();
			double mass = Mass.getAAMass(aa);
			//System.out.println("putting: " + mass);
			massFreq.put(mass, this.aaFreq.get(aa));
		}
		return massFreq;
	}
	
	
	public void printAAFreq(){
		for(Iterator<Character> it = this.aaFreq.keySet().iterator(); it.hasNext();){
			char aa = it.next();
			System.out.println(aa + "\tfreq:\t" + aaFreq.get(aa));
		}
	}
	
	public static void testGetAAFreqTable(){
		String proteinFile = "../mixture_linked/database/yeast_proteins.fasta";
		AAFreq aafreq = new AAFreq(proteinFile);
		aafreq.computeAAFreq();
		aafreq.printAAFreq();
	}
	
	public static void main(String[] args){
		testGetAAFreqTable();
	}
}
