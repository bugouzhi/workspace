package org.Spectrums;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
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
 * Estimate FDR using TDA, various info related to FDR is computed
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
	public boolean logMode = false; //some scores such as probability are usually better handle using the log-transform
	//private int scoreInd2 = 28;
	//private int rawScoreInd = 9;
	private double minScore = 0.0;
	private List<AnnotatedSpectrum> results;
	public Map<String, AnnotatedSpectrum> peptideMap;
	public Map<Integer, AnnotatedSpectrum> resultMap;
	public Map<String, AnnotatedSpectrum> proteinMap;
	private double threshold;
	public double[] thresholds = new double[]{0.01, 0.02, 0.03, 0.05, 0.1}; //we extract a set of results with standard fdr
	public int[] Prots = new int[thresholds.length];
	public int[] PSMs = new int[thresholds.length];
	public int[] peps = new int[thresholds.length];
	public boolean msplit= false;
	
	public TDAStat(String resultFile){
		this(resultFile, 7, 8, 11, 1);
	}

	public TDAStat(String resultFile, int pepInd, int protInd, int scoreInd, int sortMode){
		this.pepInd = pepInd;
		this.protInd = protInd;
		this.scoreInd = scoreInd;
		this.sortMode = sortMode;
		this.resultFile = resultFile;
		this.parseResult();
		this.getFDRByTDA();
		this.getProteinFDR();
	}
	
	
	
	public int getKeyInd() {
		return keyInd;
	}

	public void setKeyInd(int keyInd) {
		this.keyInd = keyInd;
	}

	public int getPepInd() {
		return pepInd;
	}

	public void setPepInd(int pepInd) {
		this.pepInd = pepInd;
	}

	public int getProtInd() {
		return protInd;
	}

	public void setProtInd(int protInd) {
		this.protInd = protInd;
	}

	public int getScoreInd() {
		return scoreInd;
	}

	public void setScoreInd(int scoreInd) {
		this.scoreInd = scoreInd;
	}

	public int getSortMode() {
		return sortMode;
	}

	public void setSortMode(int sortMode) {
		this.sortMode = sortMode;
	}

	public List<AnnotatedSpectrum> getResults() {
		return results;
	}
	
	public Collection<AnnotatedSpectrum> getPeptideResults() {
		return this.peptideMap.values();
	}
	
	public Collection<AnnotatedSpectrum> getProteinResults() {
		return this.proteinMap.values();
	}
	
	
	
	public int[] getPSMs() {
		return PSMs;
	}


	public void setPSMs(int[] pSMs) {
		PSMs = pSMs;
	}


	public int[] getPeps() {
		return peps;
	}


	public void setPeps(int[] peps) {
		this.peps = peps;
	}

	private boolean isSkip(String[] tokens, String line){
		return (tokens[0].startsWith("#") || tokens.length < 5 || isHeader(line));
	}
	
	private void parseResult(){
		BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(resultFile);
		//System.out.println("got reader");
		try{
			int counter = 0;
			String result = reader.readLine(); counter++;		
			results = new ArrayList<AnnotatedSpectrum>();
			resultMap = new HashMap<Integer, AnnotatedSpectrum>();
			peptideMap = new HashMap<String, AnnotatedSpectrum>();
			proteinMap = new HashMap<String, AnnotatedSpectrum>();
			while(result != null){
				String[] tokens = result.split("\\t+");
				if(isSkip(tokens, result)){
					result = reader.readLine(); counter++;
					continue;
				}
				
				
				//if(this.msplit && Double.parseDouble(tokens[11]) < 10){
				//	result = reader.readLine();
				//	continue;
				//}
				AnnotatedSpectrum s = new AnnotatedSpectrum();
				//s.scanNumber = Integer.parseInt(tokens[1]);
				s.peptide = tokens[pepInd];
				s.protein = tokens[protInd];
				s.score = Double.parseDouble(tokens[scoreInd])*sortMode;
				if(this.logMode){
					s.score = Math.log(s.score)/Math.log(10);
				}
				s.spectrumName = ""+counter;
				results.add(s);
				String key = s.peptide;// + "@";// +s.charge;
				if(peptideMap.containsKey(key)){
					AnnotatedSpectrum prev = peptideMap.get(key);
					if(s.score < prev.score){
						this.peptideMap.put(key, s);
					}
				}else{
					//if(!this.isDecoy(s.protein)){
						this.peptideMap.put(key, s);
					//}
				}
				//System.out.println("protein: " + s.protein);
				if(s.protein.contains("PROTEIN")){
					String name =  s.protein.substring(s.protein.indexOf("PROTEIN:"));
					//System.out.println("name " + name);
					name = name.split("_")[1];
					//name=name.replaceAll("_HUMAN", "");
					if(s.protein.contains("DECOY_")){
						s.protein = "DECOY_" + name;
					}else{	
						s.protein = name;
					}
				}
				//System.out.println("Protein: " + s.protein);
				key = s.protein;
				if(proteinMap.containsKey(s.protein)){
					AnnotatedSpectrum prev = proteinMap.get(s.protein);
					if(s.score < prev.score){
						prev.score = s.score;
						prev.peptide = s.peptide;
						prev.charge = s.charge;
					}
					//prev.score += s.score;
				}else{
					AnnotatedSpectrum ps = new AnnotatedSpectrum(s);
					this.proteinMap.put(key, ps);
				}
				resultMap.put(counter, s);
				result = reader.readLine(); counter++;
			}
			System.out.println("parsed result " + counter);
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
		printResultWithFDRInfo(100);
	}
	
	public void printResultWithFDRInfo(double minFDR){
		BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(resultFile);
		try{
			int counter = 0;
			String result = reader.readLine(); counter++;
			while(result != null){
				String[] tokens = result.split("\\t");
				if(isSkip(tokens, result)){
					System.out.println(result);
					result = reader.readLine(); counter++;
					continue;
				}
				if(this.resultMap.containsKey(counter)){
					AnnotatedSpectrum s = this.resultMap.get(counter);
					if((Double)s.getAnnotation().get("pepfdr") <= minFDR){
						//System.out.println(s.spectrumName + "\t" + s.peptide);
						System.out.println(result + "\t" + s.getAnnotation().get("fdr") +"\t"
								+ s.getAnnotation().get("pepfdr"));
					}
				}
				result = reader.readLine(); counter++;
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
		System.out.print("Prot:\t");
		for(int k = 0; k < this.thresholds.length; k++){
			System.out.print(this.Prots[k] + "\t");
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
	
	public void getProteinFDR(){
		List<AnnotatedSpectrum> protResults = new ArrayList();
		protResults.addAll(this.proteinMap.values());
		Collections.sort(protResults, new ScoreComparator());
		int targetCount = 0;
		int decoyCount = 0;
		int targetPep = 0;
		int decoyPep = 0;
		double protFDR = 0;
		Set<String> peps = new HashSet<String>();
		System.out.println("protein results list: " + protResults.size());
		for(int i = 0; i < protResults.size(); i++){
			AnnotatedSpectrum s = protResults.get(i);
			//System.out.println("proteins " + s.protein + "\t" + isDecoy(s.protein));
			if(isDecoy(s.protein)){
				//System.out.println("Increasing decoy count@@@");
				decoyCount++;
			}else{
				targetCount++;
			}
			protFDR = ((double)decoyCount) / targetCount;
			//System.out.println("Decoy: " + decoyCount + "\t" + "target: " + targetCount + "\t" + protFDR);
			
			for(int k = 0; k < this.thresholds.length; k++){
				if(protFDR < thresholds[k]){
					this.Prots[k] = targetCount;
				}
			}
			s.getAnnotation().put("Protfdr", protFDR);
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
	
	public static void runSwathMsplitTDA(String result, double minFDR){
		result = "../mixture_linked/UPS_Ecoli_mods_msplit_2_pep.txt";
		TDAStat stat = new TDAStat(result, 4,8,32,-1);
		stat.printResultWithFDRInfo(0.01);
		System.out.println();
	}
	
	public static void runTDA(String result, double minFDR, int pepInd, int protInd, int scoreInd, int direction){
		//result = "../mixture_linked/out.txt";
		TDAStat stat = new TDAStat(result, pepInd,protInd, scoreInd, direction);
		stat.printResultWithFDRInfo(minFDR);
		System.out.println();
	}
	
	public static void main(String[] args){
		//testTDA();
		if(args.length == 2){
			runSwathMsplitTDA(args[0], Double.parseDouble(args[1]));
		}else if(args.length == 6){
			runTDA(args[0], Double.parseDouble(args[1]), 
					Integer.parseInt(args[2]),
					Integer.parseInt(args[3]),
					Integer.parseInt(args[4]),
					Integer.parseInt(args[5]));
		}
	}
	
}
