package org.Spectrums;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
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
import UI.CommandLineParser;

/**
 * Estimate FDR using TDA, various info related to FDR is computed
 * @author Jian
 *
 */
public class TDAStat{
	private String resultFile;
	private int keyInd = 1;
	private int scanInd = -1;
	private int pepInd = 7;
	private int chargeInd = -1;
	private int protInd=8;
	private int scoreInd = 11;
	private int rtInd = -1;
	public boolean useCharge = false;
	public int sortMode = 1;  //default in natural order, but many quality scores (as oppose to stat significance) are higher the better, so use this to switch sorting mode
	public boolean logMode = false; //some scores such as probability are usually better handle using the log-transform
	//private int scoreInd2 = 28;
	//private int rawScoreInd = 9;
	private double minScore = 0.0;
	private double maxFDR = 0.01;
	private List<AnnotatedSpectrum> results;
	public Map<String, AnnotatedSpectrum> peptideMap;
	public Map<Integer, AnnotatedSpectrum> resultMap;
	public Map<String, AnnotatedSpectrum> proteinMap;
	private double threshold;
	public boolean countUniquePep=true;
	public double[] thresholds = new double[]{0.01, 0.02, 0.03, 0.05, 0.1}; //we extract a set of results with standard fdr
	public int[] Prots = new int[thresholds.length];
	public int[] PSMs = new int[thresholds.length];
	public int[] peps = new int[thresholds.length];
	public boolean msplit= false;
	public boolean hasHeader = true;
	public String fastaSeq;
	
	
	public TDAStat(String resultFile){
		this(resultFile, 7, 8, 11, 1);
	}

	public TDAStat(String resultFile, int pepInd, int protInd, int scoreInd, int sortMode){
		this(resultFile, -1, pepInd, -1, protInd, scoreInd, sortMode);
	}
	
	
	public TDAStat(String resultFile, int scanInd, int pepInd, int chargeInd, int protInd, int scoreInd, int sortMode){
		this(resultFile, scanInd, pepInd, chargeInd, protInd, scoreInd, sortMode, false);
	}
	
	/**
	 * Often it is not sure how do search engine map the peptide to protein sequence
	 * for shared peptides, which may lead to double-counting of proteins when computing protenin FDR
	 *  we try to re-map the protein at least in a consistent manner
	 */
	public TDAStat(String resultFile, int scanInd, int pepInd, int chargeInd, int protInd, int scoreInd, int sortMode, boolean useCharge){
		this(resultFile, scanInd, pepInd, chargeInd, protInd, scoreInd, sortMode, useCharge, null);
	}
	public TDAStat(String resultFile, int scanInd, int pepInd, int chargeInd, int protInd, int scoreInd, int sortMode, boolean useCharge, String fastaSeq){
		this(resultFile, scanInd, -1, pepInd, chargeInd, protInd, scoreInd, sortMode, useCharge, fastaSeq);
	}
	
	public TDAStat(String resultFile, int scanInd, int rtInd, int pepInd, int chargeInd, int protInd, int scoreInd, int sortMode, boolean useCharge, String fastaSeq){
		this(resultFile, scanInd, -1, pepInd, chargeInd, protInd, scoreInd, sortMode, useCharge, 0.01, fastaSeq);
	}
	public TDAStat(String resultFile, int scanInd, int rtInd, int pepInd, int chargeInd, int protInd, int scoreInd, int sortMode, boolean useCharge, double maxFDR, String fastaSeq){
		this.scanInd = scanInd;
		this.pepInd = pepInd;
		this.chargeInd = chargeInd;
		this.protInd = protInd;
		this.scoreInd = scoreInd;
		this.sortMode = sortMode;
		this.resultFile = resultFile;
		this.useCharge = useCharge;
		this.rtInd = rtInd;
		this.fastaSeq = fastaSeq;
		if(pepInd == 4 && protInd == 8 && scoreInd == 31){
			this.msplit = true;
		}
		this.parseResult();
		this.getFDRByTDA();
		this.maxFDR = maxFDR;
		this.getSpecCount(this.maxFDR, this.maxFDR);
		if(this.fastaSeq != null){
			this.mapProtein(this.fastaSeq);
		}
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
		if(hasHeader){
			hasHeader = false;
			return true;
		}
		return (tokens[0].startsWith("#") || tokens.length < 5 || isHeader(line) );//|| tokens[this.pepInd].contains("+"));
	}
	private String getPepKey(AnnotatedSpectrum s){
		String key;
		//System.out.println(this.useCharge);
		if(useCharge){
			key = s.peptide + "." + s.charge;
		}else{
			key = s.peptide;
		}
		return key;
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
				String[] tokens = result.split("\\t");
				if(isSkip(tokens, result)){
					result = reader.readLine(); counter++;
					continue;
				}
				
				if(this.msplit && Double.parseDouble(tokens[11]) < 10){
					result = reader.readLine(); counter++;
					continue;
				}
				//System.out.println("line is " + result);
				//System.out.println("tokens length " + tokens.length);
				AnnotatedSpectrum s = new AnnotatedSpectrum();
				if(scanInd >= 0)
					s.scanNumber = Integer.parseInt(tokens[scanInd]);
				if(this.chargeInd > 0){
					s.charge = Integer.parseInt(tokens[chargeInd]);
				}
				//s.peptide = tokens[pepInd];
				s.peptide = Utils.StringUtils.getPepSeq(tokens[pepInd]);
				s.protein = tokens[protInd];
				s.score = Double.parseDouble(tokens[scoreInd])*sortMode;
				if(this.logMode){
					s.score = Math.log(s.score)/Math.log(10);
				}
				if(this.rtInd > 0){
					s.rt = Double.parseDouble(tokens[rtInd]);
				}
				s.spectrumName = ""+counter;
				results.add(s);
				String key = getPepKey(s);
				//System.out.println("key is " + key);
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
					String name =  s.protein.substring(s.protein.indexOf("PROTEIN"));
					name = name.split("[_ ]")[1];
					//System.out.println("name " + name);
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
					//System.out.println("putting protein " + key);
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
		printResultWithFDRInfo("TDA_filtered.txt",100);
	}
	
	public void printResultWithFDRInfo(String outFile, double minPepFDR){
		printResultWithFDRInfo(outFile, minPepFDR, false);
	}
	
	//print results with FDR estimation appended at the end of each result
	public void printResultWithFDRInfo(String outFile, double minPepFDR, boolean uniquePep){
		BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(resultFile);
		try{
			int counter = 0;
			String result = reader.readLine(); counter++;
			BufferedWriter bo = new BufferedWriter(new FileWriter(outFile));
			Set<String> peps = new HashSet();
			while(result != null){
				String[] tokens = result.split("\\t");
				if(isSkip(tokens, result)){
					//System.out.println("skip");
					if(counter==1){  //write header line
						bo.write(result+"\n");
					}
					result = reader.readLine(); counter++;
					continue;
				}
				if(this.resultMap.containsKey(counter)){
					AnnotatedSpectrum s = this.resultMap.get(counter);
					if((Double)s.getAnnotation().get("pepfdr") <= minPepFDR && !s.protein.contains("DECOY_")){
						//System.out.println(s.spectrumName + "\t" + s.peptide);
						//if printing unique peptide, only print one PSM per peptide
						String pepKey = getPepKey(s);
						if(uniquePep){
							AnnotatedSpectrum pepSpect = this.peptideMap.get(pepKey);
							if(pepSpect.score == s.score && !peps.contains(getPepKey(s))){
								bo.write(result + "\t" + s.getAnnotation().get("fdr") +"\t"
									+ s.getAnnotation().get("pepfdr") + "\n");
								peps.add(pepKey);
							}
						}else{
							bo.write(result + "\t" + s.getAnnotation().get("fdr") +"fdr\t"
									+ s.getAnnotation().get("pepfdr") + "\n");
							
						}
					}
				}
				result = reader.readLine(); counter++;
			}
			bo.flush();
			bo.close();
			reader.close();
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
	
	public void getSpecCount(double minFDR, double minPepFDR){
		Map<String, Integer> specCount = new HashMap<String, Integer>();
		clearCount();
		for(int i = 0; i < results.size(); i++){
			AnnotatedSpectrum s = results.get(i);
			double fdr = (Double)s.getAnnotation().get("fdr");
			double pepFDR = (Double)s.getAnnotation().get("pepfdr");
			String key = getPepKey(s);
			if(pepFDR <= minPepFDR && fdr <= minFDR){
				if(!specCount.containsKey(key)){
					specCount.put(key, 1);
				}else{
					specCount.put(key, specCount.get(key)+1);
				}
			}
		}
		for(Iterator<String> it = this.peptideMap.keySet().iterator(); it.hasNext();){
			String pep = it.next();
			//System.out.println("pep is " + pep);
			AnnotatedSpectrum s = this.peptideMap.get(pep);
			//System.out.println("count " + specCount.get(pep));
			if(specCount.containsKey(pep)){
				s.getAnnotation().put("specCount", specCount.get(pep));
			}else{
				s.getAnnotation().put("specCount", 0);
			}
			AnnotatedSpectrum protSpect = this.proteinMap.get(s.protein);
			int count = 0;
			if(protSpect.getAnnotation().containsKey("specCount")){
				count = (Integer)protSpect.getAnnotation().get("specCount");
			}
			count+=(Integer)s.getAnnotation().get("specCount");
			protSpect.getAnnotation().put("specCount", count);
		}
		
	}
	
	private void clearCount(){
		for(Iterator<String> it = this.proteinMap.keySet().iterator(); it.hasNext();){
			String prot = it.next();
			//System.out.println("pep is " + pep);
			AnnotatedSpectrum s = this.proteinMap.get(prot);
			//System.out.println("count " + specCount.get(pep));
			s.getAnnotation().put("specCount", 0);
		}
	}
	
	public void mapProtein(String fasta){
		//System.out.println("here");
		this.fastaSeq = fasta;
		ProteinIDExtractor protID = new ProteinIDExtractor(this.results, fasta);
		Map<String, List<String>> pepMap = protID.peptideMap;
		Set<String> peps = protID.getNonSharedPeps();
		//re-mapping protein in  PSM list
		for(Iterator<AnnotatedSpectrum> it = this.resultMap.values().iterator(); it.hasNext();){
			AnnotatedSpectrum s = it.next();
			String pepKey = s.peptide;
			if(pepMap.containsKey(pepKey)){
				s.protein = pepMap.get(pepKey).get(0);
				s.getAnnotation().put("proteins", pepMap.get(pepKey));
			}else{
				s.protein = "";
			}
		}
		//re-mapping protein in peptide map
		Collection<String> newProteins = new ArrayList<String>();
		for(Iterator<AnnotatedSpectrum> it = this.peptideMap.values().iterator(); it.hasNext();){
			AnnotatedSpectrum s = it.next();
			String pepKey = s.peptide;
			if(pepMap.containsKey(pepKey)){
				s.protein = pepMap.get(pepKey).get(0);
				newProteins.add(s.protein);
				s.getAnnotation().put("proteins", pepMap.get(pepKey));
			}else{
				s.protein="";
			}
		}
		
		//re-building protein map
		Map<String, AnnotatedSpectrum> newProtMap = new HashMap<String, AnnotatedSpectrum>();
		Set<String> uniquePeps = protID.getNonSharedPeps();
		AnnotatedSpectrum dummy = new AnnotatedSpectrum();
		dummy.score = 100000;
		for(Iterator<String> it = newProteins.iterator(); it.hasNext();){
			String p = it.next();
			List<String> pepList = protID.proteinMap.get(p);
			AnnotatedSpectrum protS;
			AnnotatedSpectrum best = dummy;
			int count = 0;
			//System.out.println("processing protein " + p);
			for(Iterator<String> it2 = pepList.iterator(); it2.hasNext();){
				AnnotatedSpectrum pepS = this.peptideMap.get(it2.next());
				if(!this.countUniquePep 
						    || pepS.getAnnotation().get("proteins") == null
							|| ((List<String>)pepS.getAnnotation().get("proteins")).size() ==1){
						if(pepS.score < best.score){
							best = pepS;
						}
						count += (Integer)pepS.getAnnotation().get("specCount");
					}
			}
			if(best!=dummy){
				protS = new AnnotatedSpectrum(best);
				protS.getAnnotation().put("specCount", count);
				protS.getAnnotation().put("pepCount", protID.proteinMap.get(protS.protein).size());
				protS.getAnnotation().put("uniquePepCount", protID.getNonSharedPeps(protS.protein).size());
				//System.out.println("putting protein " + protS.protein);
				newProtMap.put(protS.protein, protS);
			}
		}
		this.proteinMap = newProtMap;
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
		return prot.contains("DECOY_") || prot.startsWith("X_") || prot.startsWith("REV_");
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
		stat.printResultWithFDRInfo(result, 0.01);
		System.out.println();
	}
	
	public static void runTDA(String result, String outFile, double minPepFDR, int pepInd, int chargeInd, int protInd, int scoreInd, int direction){
		//result = "../mixture_linked/out.txt";
		TDAStat stat = new TDAStat(result, 1, pepInd, chargeInd, protInd, scoreInd, direction);
		stat.printResultWithFDRInfo(outFile, minPepFDR);
		System.out.println();
	}
	
	
	public static void runTDAWithRT(String result, String outFile, double minPepFDR, int pepInd, int chargeInd, int protInd, int scoreInd, int direction, int RTInd1, int RTInd2){
		String firstPassOut = Utils.FileIOUtils.stripExtension(outFile) + "_withoutRT_filtered.txt";
		String resultWithRT = Utils.FileIOUtils.stripExtension(result)+ "_withRT_allresult.txt";
		TDAStat stat = new TDAStat(result, 1, pepInd, chargeInd, protInd, scoreInd, direction);
		stat.printResultWithFDRInfo(firstPassOut, minPepFDR);
		//RTAligner.runRTAlign(firstPassOut, result, resultWithRT, RTInd1, RTInd2, pepInd, chargeInd);
		System.out.println();
	}
	
	//test running TDA on a dir of results file
	public static void runTDADir(String resultDir, String outFile, double minPepFDR, int pepInd, int chargeInd, int protInd, int scoreInd, int direction){
		String combined = outFile+".combined.tmp";
		Utils.FileIOUtils.mergeResultWithHeader(resultDir, combined);
		TDAStat stat = new TDAStat(combined, 1, pepInd, chargeInd, protInd, scoreInd, direction, true);
		stat.printResultWithFDRInfo(outFile, minPepFDR, true);
		System.out.println();
		File tmp = new File(combined);
		try{		
			Files.deleteIfExists(tmp.toPath());
		}catch(Exception e){
			System.out.println(e.getMessage());
		}
	}
	
	
	public static void main(String[] args){
		//testTDA();
		//args[0] = "SWATH-MSPLIT";
		CommandLineParser parser = new CommandLineParser(args);
		if(args.length == 2){
			runSwathMsplitTDA(args[0], Double.parseDouble(args[1]));
		}else if(args.length == 7){
			runTDA(args[0], args[1], Double.parseDouble(args[2]), 
					1,
					Integer.parseInt(args[3]),
					Integer.parseInt(args[4]),
					Integer.parseInt(args[5]),
					Integer.parseInt(args[6]));
		}else if(args.length == 4){
			if(args[0].equals("SWATH-MSPLIT")){
				//args[1] = "../mixture_linked/SWATH/testRTConstraints/Human_500nglysate_QTOF5600_tppSpstconsensus_lib_swathmsplit_out_rtFiltered.txt";
				//args[2] = "../mixture_linked/SWATH/testRTConstraints/filtered.txt";
				//args[3] = "0.01";
				runTDA(args[1], args[2], Double.parseDouble(args[3]), 4, 6, 8, 32 , -1);
			}
			if(args[0].equals("MSGFDB-Dir")){
				//args[1] = "../mixture_linked/SWATH-MSPLITv1.0/test_DDA_results/";
				//args[2] = "../mixture_linked/SWATH-MSPLITv1.0/test_DDA_combinedfiltered.txt";
				//args[3] = "0.01";
				runTDADir(args[1], args[2], Double.parseDouble(args[3]), 7, 6, 8, 12 , 1);
			}else if(args[0].equals("Generic-Dir")){
				runTDADir(args[1], args[2], Double.parseDouble(args[3]),
						parser.getInteger(4), parser.getInteger(5), parser.getInteger(6), parser.getInteger(7), parser.getInteger(8));
			}
		}
	}
	
}
