package org.Spectrums;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;

/**
 * Compute various TDA related info for a set of PSMs
 * @author Jian
 *
 */
public class TDAStat{
	private String resultFile;
	private int keyInd = 1;
	private int pepInd = 7;
	private int protInd=8;
	private int scoreInd = 11;
	public int sortMode = 1;  //default in natural order, but many quality scores (as oppose to stat significance) are higher the better, so use this to switch sorting mode
	//private int scoreInd2 = 28;
	//private int rawScoreInd = 9;
	private double minScore = 0.0;
	private List<AnnotatedSpectrum> results;
	private Map<Integer, AnnotatedSpectrum> resultMap;
	private double threshold;
	public double[] thresholds = new double[]{0.01, 0.02, 0.03, 0.05, 0.1};
	private int[] PSMs = new int[thresholds.length];
	private int[] peps = new int[thresholds.length];
	
	public TDAStat(String resultFile){
		this.resultFile = resultFile;
		this.parseResult();
		this.getFDRByTDA();
	}
	
	
	private void parseResult(){
		BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(resultFile);
		//System.out.println("got reader");
		try{
			String result = reader.readLine();
			int counter = 0;
			results = new ArrayList<AnnotatedSpectrum>();
			resultMap = new HashMap<Integer, AnnotatedSpectrum>();
			while(result != null){
				if(isHeader(result)){
					result = reader.readLine();
					continue;
				}
				String[] tokens = result.split("\\t");
				AnnotatedSpectrum s = new AnnotatedSpectrum();
				s.scanNumber = Integer.parseInt(tokens[1]);
				s.peptide = tokens[pepInd];
				s.protein = tokens[protInd];
				s.score = Double.parseDouble(tokens[scoreInd]);
				s.spectrumName = ""+counter;
				results.add(s);
				resultMap.put(counter, s);
				counter++;
				result = reader.readLine();
			}
			reader.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	/**
	 * Print out the original results with FDR info added
	 */
	public void printResultWithFDRInfo(){
		BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(resultFile);
		try{
			String result = reader.readLine();
			int counter = 0;
			while(result != null){
				if(isHeader(result)){
					System.out.println(result);
					result = reader.readLine();
					continue;
				}
				AnnotatedSpectrum s = this.resultMap.get(counter);
				System.out.println(result + "\t" + s.getAnnotation().get("fdr") +"\t"
						+ s.getAnnotation().get("pepfdr"));
				result = reader.readLine();
				counter++;
			}
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public void getSummary(){
		System.out.println("TDA-Result-Summary:");
		System.out.print(" FDR:\t");
		for(int k = 0; k < this.thresholds.length; k++){
			System.out.print(this.thresholds[k] + "\t");
		}
		System.out.println();
		System.out.print("PSMs:\t");
		for(int k = 0; k < this.thresholds.length; k++){
			System.out.print(this.PSMs[k] + "\t");
		}

		System.out.println();
		System.out.print("Peps:\t");
		for(int k = 0; k < this.thresholds.length; k++){
			System.out.print(this.peps[k] + "\t");
		}
		System.out.println();

		
	}
	
	public void getFDRByTDA(){
		Collections.sort(results, new ScoreComparator());
		int targetCount = 0;
		int decoyCount = 0;
		int targetPep = 0;
		int decoyPep = 0;
		double FDR = 0;
		double pepFDR = 0;
		Set<String> peps = new HashSet<String>();
		for(int i = 0; i < results.size(); i++){
			AnnotatedSpectrum s = results.get(i);
			//System.out.println("proteins " + s.protein + "\t" + isDecoy(s.protein));
			if(isDecoy(s.protein)){
				decoyCount++;
				if(!peps.contains(s.peptide)){
					peps.add(s.peptide);
					decoyPep++;
				}
			}else{
				targetCount++;
				if(!peps.contains(s.peptide)){
					peps.add(s.peptide);
					targetPep++;
				}
			}

		
			FDR = ((double)decoyCount) / targetCount;
			pepFDR =  ((double)decoyPep) / targetPep;
			
			for(int k = 0; k < this.thresholds.length; k++){
				if(FDR < thresholds[k]){
					this.PSMs[k] = targetCount;
				}
			}
			
			for(int k = 0; k < this.thresholds.length; k++){
				if(pepFDR < thresholds[k]){
					this.peps[k] = targetPep;
				}
			}

			s.getAnnotation().put("fdr", FDR);
			s.getAnnotation().put("pepfdr", pepFDR);
		}
	}
		
	protected double getThreshold(List<Double> target, List<Double> decoy, double fdr){
		int totalTarget = target.size();
		int totalDecoy = decoy.size();
		System.out.println("Total numbers of targets: " + target.size());
		System.out.println("Total numbers of decoys: " + decoy.size());
		Collections.sort(target);
		Collections.sort(decoy);
		int i = 0, j = 0; 
		while(i < target.size() && j < decoy.size()){
			//System.out.println("target: " + (target.size()-i) +"\t" + "decoy: " + (decoy.size()-j) + "\t" + target.get(i) + "\t" + decoy.get(j));
			
			if(target.get(i) <= decoy.get(j)){
				i++;
			}else{
				j++;
			}
			if((double)(totalDecoy - j ) / (double)(totalTarget - i ) < fdr){
				break;
			}
		}
		System.out.println("threshold is: " + target.get(i-1) + " Accepted targets: " + (totalTarget - i));
		if(i == target.size()){
			i=i-1;
		}
		return target.get(i);
	}
	
	
	protected boolean isDecoy(String prot){
		return prot.contains("DECOY_") || prot.startsWith("X_");
	}
	
	protected boolean isHeader(String result){
		return result.startsWith("#");
	}
	
	class ScoreComparator implements Comparator{
		public ScoreComparator(){
			
		}
		public int compare(Object s1, Object s2) {
			return (int)compare((Spectrum)s1, (Spectrum)s2);
		}
		public double compare(Spectrum s1, Spectrum s2) {
			if(s1.score < s2.score){
				return -1;
			}else if(s1.score == s2.score){
				return 0;
			}else{
				return 1;
			}
		}		
	}
	
	public static void testTDA(){
		String resultFile = "../mixture_linked/MSGFDBSearches/14341_UPS1-400fm_IDA_msgfdb.txt";
		TDAStat stat = new TDAStat(resultFile);
		stat.printResultWithFDRInfo();
		System.out.println();
		stat.getSummary();
	}
	
	public static void main(String[] args){
		testTDA();
	}
	
}
