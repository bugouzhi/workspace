package mixgf;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.Spectrums.AAFreq;
import org.Spectrums.LargeSpectrumLibIterator;
import org.Spectrums.Mass;
import org.Spectrums.MixturePeakScoreLearner;
import org.Spectrums.PeakComparator;
import org.Spectrums.Peptide;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.SimpleProbabilisticScorer;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumComparator;
import org.Spectrums.SpectrumLib;
import org.Spectrums.TheoreticalSpectrum;

import IO.MZXMLReader;
import Utils.ArrayUtils;
import Utils.StringUtils;

/**
 * Compute spectral probability for mixture spectrum using
 * the generating function approach
 * @author Jian
 *
 */

public class MXGF {
	public static double MAXMASS = 3000;
	public static double MAXSCORE = 250;
	public static double MINSCORE = -200;
	Map<Double, Double> massFreqTable;
	int[] validAAMass;
	double[] aaFreq; //note aa freq table should be matched to the validAAMass;
	double maxMass = 0;
	double maxMass2 = 0;
	double minScore = 0;
	double minScore2 = 0;
	double maxScore =0;
	double maxScore2 = 0;
	double[] scoreHi;
	double[] scoreLow;
	double[][][] table; //DP table
	double[] pepCounts;
	double resolution = 1.0;
	int maxAAMass= 187;
	int currentIndex;
	int currentIndex2;
	boolean DEBUG = false;
	double scaleFactor = 0.9995;
	double tolerance = 0.5;
	
	public MXGF(){
		initialize();
	}
	
	public MXGF(double maxMass, double maxScore){
		this(maxMass,maxMass, maxScore);
	}
	
	public MXGF(double maxMass, double maxMass2, double maxScore){
		this(maxMass, maxMass2, maxScore, null);
	}
	
	public MXGF(double maxMass, double maxScore, Map<Double, Double> aaFreq){
		this(maxMass, maxMass, maxScore, aaFreq);
	}
	
	public MXGF(double maxMass, double maxMass2, double maxScore, Map<Double, Double> aaFreq){
		this(maxMass, maxMass2, maxScore, aaFreq, 1.0);
	}
	
	public MXGF(double maxMass, double maxMass2, double maxScore, Map<Double, Double> aaFreq, double resolution){
		this.resolution = resolution;
		this.maxMass = Math.ceil(maxMass/this.resolution);
		this.maxMass2 = Math.ceil(maxMass2/this.resolution);
		this.maxScore = maxScore;
		this.massFreqTable = aaFreq;
		initialize();
	}
	
	private void initialize(){
		if(this.maxMass == 0)
			maxMass = MAXMASS/this.resolution;
		if(this.maxMass2 == 0)
			maxMass2 = MAXMASS/this.resolution;
		if(this.maxScore == 0)
			maxScore = MAXSCORE;
		if(this.minScore == 0)
			minScore = MINSCORE;
		this.scoreHi = new double[((int)Math.ceil(maxMass/this.resolution))+2];
		this.scoreLow = new double[((int)Math.ceil(maxMass2/this.resolution))+2];
		this.getAAStat();
		double max = maxMass > maxMass2 ? maxMass : maxMass2;
		this.table = new double[(int)Math.ceil(max)][2][(int)Math.ceil(maxScore-minScore+1)];
		if(DEBUG){
			System.out.println("scoreHi size: " + this.scoreHi.length);
			System.out.println("scoreLo size: " + this.scoreLow.length);
			System.out.println("table size: " + this.table.length + "\t" + this.table[0].length);
		}
	}
	
	private void getAAStat(){
		int validAA = 0;
		for(int i = 0; i < Mass.aaMap.length; i++){
			if(Mass.aaMap[i] > 0 && Mass.aaMap[i] != 111.0316){
				validAA++;
			}
		}
		this.validAAMass = new int[validAA];
		this.aaFreq = new double[validAA];
		int j = 0;
		for(int i = 0; i < Mass.aaMap.length; i++){
			if(Mass.aaMap[i] > 0 && Mass.aaMap[i] != 111.0316){
				if(DEBUG){
					System.out.println("valid aa " + Math.round(Mass.aaMap[i]*scaleFactor/this.resolution));
				}
				this.validAAMass[j] = (int)Math.round(Mass.aaMap[i]*scaleFactor/this.resolution);
				if(this.massFreqTable != null){
					System.out.println("mass : " + Mass.aaMap[i]);
					if(massFreqTable.containsKey(Mass.aaMap[i])){
						this.aaFreq[j++] = this.massFreqTable.get(Mass.aaMap[i]);
					}else{
						this.aaFreq[j++] = 0.0; //does it make sense to have zero here or should we use very small number instead?
					}
				}else{
					this.aaFreq[j++] = 0.05; //assume uniform model
				}
			}
		}		
	}
	

	//compute probability for one peptide only i.e. compute table(m1,0,s)
	public void initTable(){
		//initialize table
		if(DEBUG) System.out.println("score length: " + this.scoreHi.length + "\t" + this.scoreLow.length);
		double maxScore1 = Math.ceil(this.computeMaxScore(this.scoreHi))+5;
		double minScore1 = Math.floor(this.computeMinScore(this.scoreHi))-5;
		double maxScore2 = Math.ceil(this.computeMaxScore(this.scoreLow))+5;
		double minScore2 = Math.floor(this.computeMinScore(this.scoreLow))-5;
		this.maxScore = maxScore1 > maxScore2 ? maxScore1 : maxScore2;
		this.minScore = minScore1 < minScore2 ? minScore1 : minScore2;
		//System.out.println("max Score is: " + this.maxScore);
		//System.out.println("min Score is: " + this.minScore);
		double max = maxMass > maxMass2 ? maxMass : maxMass2;
		if(DEBUG) System.out.println("table size: " + (int)Math.ceil(max) + "\t" + (int)Math.ceil(maxScore-minScore+1));
		this.table = new double[(int)Math.ceil(max)][2][(int)Math.ceil(maxScore-minScore+1)];		
	}
	
	public void computeHi(){
		int offSet = (int)(Math.round(-minScore));
		table[0][0][offSet] = 1;
		int minS = (int)Math.floor(minScore);
		//System.out.println("offset: " + scoreOffSet);
		//System.out.println("maxmass " + maxMass);
		for(int m1 = 0; m1 < maxMass; m1++){
			for(int s = minS; s < maxScore; s++){
				for(int a = 0; a < this.validAAMass.length; a++){
					int prevMass = m1 - this.validAAMass[a];
					int prevScore = (int)(s - this.scoreHi[m1]);  //check rounding errors 
 					double prevCount = 0;
					if(prevMass >= 0 && prevScore > minScore && prevScore < maxScore){
						prevScore = prevScore - minS;
						prevCount = this.table[prevMass][0][prevScore];
						table[m1][0][s - minS] += prevCount*0.05;//*this.aaFreq[a];
						if(DEBUG) {
							if(table[m1][0][s - minS] > 0){
								//System.out.println("prob at: " + m1 + "\t" + table[m1][0][s-minS]);
							}
						}
					}
				}
			}
		}
		//this.computeNumPeptides();
	}
	
	
	/**
	 * Compute the conditional probabilities
	 * @param masses
	 */
	public void computeCondPair(int[] masses){
		this.initTable();
		int offSet = (int)(Math.round(-minScore));
		table[0][1][offSet] = 1;
		int minS = (int)Math.floor(minScore);
		double[] oldScoresHi = new double[masses.length];
		double[] oldScoresLow = new double[masses.length];
		for(int i = 0; i < masses.length; i++){
			if(masses[i] <= this.maxMass2){
				oldScoresLow[i] = this.scoreLow[masses[i]];
				if(masses[i] <= this.maxMass && masses[i] <= this.maxMass2){
					if(this.scoreHi[masses[i]] >= this.scoreLow[masses[i]]){
						//this.scoreLow[masses[i]]  = 0.0;
					}
				}
			}
			if(masses[i] <= this.maxMass && masses[i] <= this.maxMass2){
				oldScoresHi[i] = this.scoreHi[masses[i]];
				if(this.scoreHi[masses[i]] < this.scoreLow[masses[i]]){
					//this.scoreHi[masses[i]]  = 0.0;
				}
			}
		}
		this.computeHi();
		//System.out.println("offset: " + scoreOffSet);
		for(int m1 = 0; m1 < maxMass2; m1++){
			for(int s = minS; s < maxScore; s++){
				for(int a = 0; a < this.validAAMass.length; a++){
					int prevMass = m1 - this.validAAMass[a];
					int prevScore = (int)(s - this.scoreLow[m1]);
 					double prevCount = 0;
					if(prevMass >= 0 && prevScore > minScore && prevScore < maxScore){
						prevScore = prevScore - minS;
						prevCount = this.table[prevMass][1][prevScore];
						table[m1][1][s - minS] += prevCount*0.05;//*this.aaFreq[a];
					}
				}
			}
		}
		for(int i = 0; i < masses.length; i++){
			if(masses[i] <= this.maxMass2){
				this.scoreLow[masses[i]] = oldScoresLow[i];
			}
			
			if(masses[i] <= this.maxMass){
				this.scoreHi[masses[i]] = oldScoresHi[i];
			}
		}
	}
	
	/**
	 * Compute joint probability using a dynammic programming approach. Note
	 * since it is too big to fit the full 3D DP matrix in memory we use
	 * a "sliding-window" approach, so for computing the DP at a particular location
	 * we only keep the DP matrix before that is require to compute the current location
	 * (i.e. the largest amino acid mass)
	 */
	public void computePair(){
		this.table = new double[this.maxAAMass+1][this.maxAAMass+1][(int)Math.ceil(maxScore-minScore+1)];
		//System.out.println("max Score is: " + this.maxScore);
		//System.out.println("min Score is: " + this.minScore);
		int offSet = (int)(-minScore);
		table[0][0][offSet] = 1;  
		int minS = (int)Math.floor(minScore);
		//System.out.println("offset: " + scoreOffSet);
		currentIndex=-1*this.maxAAMass; //don't need to loop around until reach first 187
		currentIndex2=-1*this.maxAAMass;
		int m1=0, m2 = 0, currentMass=1, currentMass2=1;
		double aaProb = 0.05;
		double max = maxMass > maxMass2 ? maxMass : maxMass2;
		//System.out.println(maxMass + "\t" + maxMass2 + "\t" + this.scoreHi.length + "\t" + this.scoreLow.length);
		while(currentMass <= max && currentMass2 <= max){
				//System.out.println("current mass: " + currentMass + "\t" + currentMass2 + "\t" + currentIndex + "\t" + currentIndex2);
				//case I: if m2 < m1
				m1 = currentMass;
				m2 = currentMass2 - this.maxAAMass;
				m2 = m2 > 0 ? m2 : 0;
				if(currentMass <= maxMass){
				for(; m2 < m1 && m2 < maxMass2; m2++){
					for(int s = minS; s < maxScore; s++){
						table[getIndex(currentMass, 0, currentIndex)][getIndex(currentMass2, currentMass2-m2, currentIndex2)][s - minS]=0;
						for(int a = 0; a < this.validAAMass.length; a++){
							int prevMass = m1 - this.validAAMass[a];
							int prevScore = (int)(s - this.scoreHi[m1]);
							double prevCount = 0;
							int count = 0;
							if(prevMass >= 0 && prevScore > minScore && prevScore < maxScore){
								prevScore = prevScore - minS;
								prevCount = this.table[getIndex(currentMass, this.validAAMass[a], currentIndex)][getIndex(currentMass2, currentMass2-m2, currentIndex2)][prevScore];
								table[getIndex(currentMass, 0, currentIndex)][getIndex(currentMass2, currentMass2-m2, currentIndex2)][s - minS] += prevCount*aaProb;//this.aaFreq[a];
								count++;
							}
						}
					}
				}
				}
				//case II: if m1 < m2
				m2 = currentMass2;
				m1 = currentMass - this.maxAAMass;
				m1 = m1 > 0 ? m1 : 0;
				if(currentMass2 < maxMass2){
				for(; m1 < m2 && m1 < maxMass; m1++){
					for(int s = minS; s < maxScore; s++){
						table[getIndex(currentMass, currentMass-m1, currentIndex)][getIndex(currentMass2, 0, currentIndex2)][s - minS] = 0;
						for(int a = 0; a < this.validAAMass.length; a++){
							int prevMass = m2 - this.validAAMass[a];
							int prevScore = (int)(s - this.scoreLow[m2]);
							double prevCount = 0;
							int count = 0;
							if(prevMass >= 0 && prevScore > minScore && prevScore < maxScore){
								prevScore = prevScore - minS;
								prevCount = this.table[getIndex(currentMass, currentMass-m1, currentIndex)][getIndex(currentMass2, this.validAAMass[a], currentIndex2)][prevScore];
								table[getIndex(currentMass, currentMass-m1, currentIndex)][getIndex(m2, 0, currentIndex2)][s - minS] += prevCount*aaProb;//this.aaFreq[a];
								count++;
							}
						}
					}
				}
				}
				//case III: if m1 == m2
				if(currentMass == currentMass2){
				m1 = currentMass; m2 = currentMass2;
				for(int s = minS; s < maxScore; s++){
					table[getIndex(currentMass, 0, currentIndex)][getIndex(currentMass2, 0, currentIndex2)][s - minS]=0;
					for(int a = 0; a < this.validAAMass.length; a++){
						for(int a2 = 0; a2 < this.validAAMass.length; a2++){
							int prevMass1 = m1 - this.validAAMass[a];
							int prevMass2 = m2 - this.validAAMass[a2];
							int prevScore1 = (int)(s - this.scoreHi[m1]);
							int prevScore2 = (int)(s - this.scoreLow[m2]);
							//int prevScore = prevScore1 > prevScore2 ? prevScore2 : prevScore1;
							int prevScore = (int)(s - this.scoreHi[m1]-this.scoreLow[m2]);
							double prevCount = 0;
							if(prevMass1 >= 0 && prevMass2 >=0 && prevScore > minScore && prevScore < maxScore){
								prevScore = prevScore - minS;
								prevCount = this.table[getIndex(currentMass, this.validAAMass[a], currentIndex)][getIndex(currentMass2, this.validAAMass[a2], currentIndex2)][prevScore];
								table[getIndex(currentMass, 0, currentIndex)][getIndex(currentMass2, 0, currentIndex2)][s - minS] += prevCount*aaProb*aaProb;//this.aaFreq[a];
							}
						}
					}
				}
				}
			if(currentMass <= maxMass){	
				currentIndex=nextCurrentIndex(currentIndex);
				currentMass++;
			}
			if(currentMass2 <= maxMass2){
				currentIndex2=nextCurrentIndex(currentIndex2);
				currentMass2++;
			}			
		}
		
		//last loop don't need to advance currentIndex, so can query the table correctly
		if(currentIndex == 0){   
			currentIndex = this.table.length-1;
		}else{
			currentIndex--;
		}
		if(currentIndex2 == 0){   
			currentIndex2 = this.table.length-1;
		}else{
			currentIndex2--;
		}
	}
	
	private int getIndex(int mass, int jump, int currentIndex){
		int index = currentIndex - jump;
		if(currentIndex < 0){
			index = mass-jump;
		}else if(index < 0){
			index = this.maxAAMass + index +1;
		}
		return index;
	}
	
	private int nextCurrentIndex(int currentIndex){
		 if(currentIndex == this.table.length-1){
			return 0;
		}else{
			return currentIndex+1;
		}
	}
	
	
	
	public void computeNumPeptides(){
		this.pepCounts = new double[table.length];
		//System.out.println("offset: " + scoreOffSet);
		pepCounts[0]=1;
		for(int m1 = 0; m1 < maxMass; m1++){
			for(int a = 0; a < this.validAAMass.length; a++){
				int prevMass = m1 - this.validAAMass[a];
				double prevCount = 0;
				if(prevMass >= 0){
					prevCount = pepCounts[prevMass];//*this.aaFreq[a];
					pepCounts[m1]+=prevCount;
				}
			}
		}
	}
	
	
	/**
	 * Compute the maximum score that can be achieved using the current scoring model
	 * @param scores
	 * @return
	 */
	public double computeMaxScore(double[] scores){
		double[] maxScore = new double[scores.length];
		//System.out.println("offset: " + scoreOffSet);
		maxScore[0]=0;
		for(int m1 = 1; m1 < scores.length; m1++){
			double max = -1000;
			for(int a = 0; a < this.validAAMass.length; a++){
				int prevMass = m1 - this.validAAMass[a];
				if(prevMass >= 0){
					 double current = maxScore[prevMass];//*this.aaFreq[a];
					 max = current > max ? current : max;
				}
			}
			maxScore[m1]+=max+scores[m1];
		}
		
		double max = -1000;
		for(int i = 0; i < maxScore.length; i ++){
			max = maxScore[i] > max ? maxScore[i] : max;
			//System.out.println("max at : " + i +"\t" + maxScore[i]);
		}
		if(DEBUG) System.out.println("max score at mass: " + maxMass + " is: " + maxScore[maxScore.length-1] + " the max is: " + max);
		return max;
	}
	
	/**
	 * Compute the minimum score that can be achieved using the current scoring model
	 * @param scores
	 * @return
	 */
	public double computeMinScore(double scores[]){
		double[] minScore = new double[scores.length];
		//System.out.println("offset: " + scoreOffSet);
		minScore[0]=0;
		for(int m1 = 1; m1 < scores.length; m1++){
			double min = 1000;
			for(int a = 0; a < this.validAAMass.length; a++){
				int prevMass = m1 - this.validAAMass[a];
				if(prevMass >= 0){
					 double current = minScore[prevMass];//*this.aaFreq[a];
					 min = current < min ? current : min;
				}
			}
			minScore[m1]+=min+scores[m1];
		}
		double min = 1000;
		for(int i = 0; i < minScore.length; i ++){
			min = minScore[i] < min ? minScore[i] : min;
		}
		return min;
	}
	
	private int getScoreIndex(int score){
		return 0;
	}
	
	//for single peptide-hi
	public double getSpecProb(double score, double mass){
		return getSpecProb(score, mass, 0);
	}
	
	/**
	 * Get the spectral probability at a particular mass and score position
	 * 
	 * @param score1
	 * @param score2
	 * @param mass1
	 * @param mass2
	 * @return
	 */
	public double[] getSpecProb(double score1, double score2, double mass1, double mass2){
		return new double[]{getSpecProb(score1, mass1, 0), 
				getSpecProb(score2, mass2, 1)};
	}
	
	public double getJointProb(double score, double mass1, double mass2){
		int s = (int)Math.round(score-minScore);
		int m1 = (int)Math.round(this.scaleMass(mass1)/this.resolution);
		int m2 = (int)Math.round(this.scaleMass(mass2)/this.resolution);
		System.out.println("massses for prob: " + m1 + "\t" + m2);
		if(this.maxMass - mass1 > this.maxAAMass || this.maxMass2 - mass2 > this.maxAAMass){
			System.err.println("warning: querying mass out of range: " + mass1 + "\t" + mass2);
			return 0.0;
		}
		m1 = getIndex((int)maxMass, (int)(maxMass-m1), currentIndex);
		m2 = getIndex((int)maxMass2, (int)(maxMass2-m2), currentIndex2);
		double total = 0;
		double higher = 0.0;
		double highest = 0.0;
		for(int currentS = 0; currentS < maxScore-minScore; currentS++){
			if(currentS >= s){
				higher+=table[m1][m2][currentS];
			}
			total += table[m1][m2][currentS];
			if(table[m1][m2][currentS] > 0){
				highest=table[m1][m2][currentS];
			}
			//System.out.println((currentS + minScore) + "\t" + table[m1][m2][currentS]);
		}
		if(higher == 0){
			higher=highest; //if higher is zero we use highest
		}
		//System.out.println("total possible candidates: " + this.pepCounts[m1]);
		//System.out.println("total: " + total + "\t" + "higher " + higher);
		return higher;///total;
	}
	
	private double getSpecProb(double score, double mass, int index){
		int s = (int)Math.round(score-minScore);
		int m1 = (int)Math.round(this.scaleMass(mass)/this.resolution);
		//System.out.println("m1 " + m1 + "\t" + score + "\t" + s);
		double total = 0;
		double higher = 0.0;
		double highest = 0.0;
		for(int currentS = 0; currentS < maxScore-minScore; currentS++){
			if(currentS >= s){
				higher+=table[m1][index][currentS];
			}
			total += table[m1][index][currentS];
			if(table[m1][index][currentS] > 0){
				highest=table[m1][index][currentS];
			}
			//System.out.println((currentS + minScore) + "\t" + table[m1][index][currentS]);
		}
		if(higher == 0){
			higher=highest; //if higher is zero we use highest
		}
		//System.out.println("total possible candidates: " + this.pepCounts[m1]);
		//System.out.println("total: " + total + "\t" + "higher " + higher);
		return higher;///total;
	}
		
	public double scaleMass(double mass){
		return scaleFactor*mass;
	}
	public double getMaxMass() {
		return maxMass;
	}

	public void setMaxMass(double maxMass) {
		this.maxMass = maxMass;
	}

	public double getMaxScore() {
		return maxScore;
	}

	public void setMaxScore(double maxScore) {
		this.maxScore = maxScore;
	}

	public double[] getScoreHi() {
		return scoreHi;
	}

	public void setScoreHi(double[] scoreHi) {
		for(int i = 0; i < scoreHi.length; i++){
			//System.out.println(i + "\tscore: " + scoreHi[i]);
			if(i < this.scoreHi.length)
				this.scoreHi[i] = scoreHi[i];
		}
	}

	public double[] getScoreLow() {
		return scoreLow;
	}

	public void setScoreLow(double[] scoreLow) {
		for(int i = 0; i < scoreLow.length; i++){
			//if(DEBUG) System.out.println(i + "\tscore: " + scoreLow[i]);
			if(i < this.scoreLow.length)
				this.scoreLow[i] = scoreLow[i];
		}
	}

	public double[][][] getTable() {
		return table;
	}

	public void setTable(double[][][] table) {
		this.table = table;
	}
	
	protected static void testComputeMSGF(){
		String spectrumFile = "../mixture_linked/yeast_data/klc_010908p_yeast-digest.mzXML";
		String scorerFile = "../mixture_linked/yeast_single_peptide_model.o";
		String proteinsFile = "../mixture_linked/database/yeast_proteins.fasta";
		String resultFile = "../mixture_linked/testAnnotation.txt";
		String prmFile = "../mixture_linked/MSGFSearches/testPRMOut_noenzy.txt";
		List<String> resultLines = Utils.FileIOUtils.createListFromFile(resultFile);
		AAFreq aafreq = new AAFreq(proteinsFile);
		MZXMLReader iterator = new MZXMLReader(spectrumFile);
		SpectrumLib lib = new SpectrumLib(prmFile, "MGF");
		PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		((SimpleProbabilisticScorer)scorer).setMinMatchedPeak(0);
		List<Spectrum> prmList = lib.getAllSpectrums();
		for(Iterator<String> it = resultLines.iterator(); it.hasNext();){
			String line = it.next();
			String[] tokens = line.split("\\s+");
			Spectrum s = iterator.getSpectrum(Integer.parseInt(tokens[3]));
			if(s.scanNumber != 1004){
				//continue;
			}
			s.windowFilterPeaks(5, 25);
			s.computePeakRank();
			String peptide = tokens[7];
			if(peptide.startsWith("r")){
				peptide = peptide.substring(1);
			}
			Peptide p = new Peptide(peptide);
			System.out.println("peptide is: " + p + "\t" + p.getParentmass() +"\t" + p.getCharge() + "\tspectrum: " + s.parentMass + "\t" + s.charge);
			//s.parentMass = p.getParentmass();
			//s.charge = p.getCharge();
			TheoreticalSpectrum t = new TheoreticalSpectrum(p);
			double score = scorer.compare(t, s);
			double[][] base = t.computeBaseMass(p.getPeptide(), p.getPos(), p.getPtmmasses());
			double[] scores = null;
			for(int i = 0; i < prmList.size(); i++){
				Spectrum prm = prmList.get(i);
				if(prm.scanNumber == Integer.parseInt(tokens[3])){
					scores = new double[prm.getPeak().size()+1];
					for(int j = 1; j < prm.getPeak().size(); j++){
						scores[j] = prm.getPeak().get(j-1).getIntensity();
					}
				}
			}
			double totalScore = 0;
			base[0][base[0].length-1]=p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER;
			for(int i = 0; i < base[0].length; i++){
				int massIndex = (int)Math.round((0.9995*base[0][i]));
				double curr=0;
				if(massIndex > 0 && massIndex < scores.length)
					curr= scores[massIndex];
				System.out.println("scored: " + massIndex +"\t" + curr);
				totalScore += curr;
			}
			MXGF mixgf = new MXGF(t.parentMass*t.charge, 200, aafreq.getMassFreq());
			mixgf.setScoreHi(scores);
			mixgf.computeHi();
			double prob = mixgf.getSpecProb(totalScore, (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
			System.out.println("score is: " + score + "\tPRM total score:\t" + totalScore);
			System.out.println(line + "\t" + score +"\t" + totalScore + "\t" + prob);
			//return;
		}
		
	}
	
	protected static void testComputeSingleMXGF(){
		String spectrumFile = "../mixture_linked/yeast_data/klc_010908p_yeast-digest.mzXML";
		String scorerFile = "../mixture_linked/yeast_NIST_single_peptide_model.o";
		//String scorerFile = "../mixture_linked/mixtures_alpha_generic.o";
		String proteinsFile = "../mixture_linked/database/yeast_proteins.fasta";
		String resultFile = "../mixture_linked/testAnnotation2.txt";
		List<String> resultLines = Utils.FileIOUtils.createListFromFile(resultFile);
		AAFreq aafreq = new AAFreq(proteinsFile);
		MZXMLReader iterator = new MZXMLReader(spectrumFile);	
		PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		//PeakComparator comp = MixturePeakScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		((SimpleProbabilisticScorer)scorer).setMinMatchedPeak(0);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int count = 0;
		for(Iterator<String> it = resultLines.iterator(); it.hasNext();){
			String line = it.next();
			String[] tokens = line.split("\\s+");
			Spectrum s = iterator.getSpectrum(Integer.parseInt(tokens[2]));
			if(s.scanNumber != 593){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			String peptide = tokens[10];
			if(peptide.startsWith("r")){
				peptide = peptide.substring(1);	
			}
			String peptide2 = tokens[12];
			if(peptide.startsWith("r")){
				peptide = peptide.substring(1);	
			}
			
			Peptide p = new Peptide(peptide);
			Peptide p2 = new Peptide(peptide2);
			System.out.println("peptide is: " + p + "\t" + p.getParentmass() +"\t" + p.getCharge() + "\t" + p2 + "\t" + p2.getParentmass() + "\t" + p2.getCharge()+ "\tspectrum: " + s.parentMass + "\t" + s.charge);
			s.parentMass = p.getParentmass();
			s.charge = p.getCharge();
			TheoreticalSpectrum t = new TheoreticalSpectrum(p);
			//((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).mode = 1;
			PRMSpectrum prmSpect = new PRMSpectrum(s, p.getCharge(), scorer, 1.0);
			double score = scorer.compare(t, s);
			s.parentMass = p2.getParentmass();
			s.charge = p2.getCharge();
			//((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).mode = 2;
			TheoreticalSpectrum t2 = new TheoreticalSpectrum(p2);
			PRMSpectrum prmSpect2 = new PRMSpectrum(s, p2.getCharge(), scorer, 1.0);
			double score2 = scorer.compare(t2, s);
			double[] scores = prmSpect.getScoredSpectrum(p.getParentmass()*p.getCharge());
			double[] scores2 = prmSpect2.getScoredSpectrum(p2.getParentmass()*p2.getCharge());
			PRMSpectrumComparator sComp = new PRMSpectrumComparator();
			double[][] base = prmSpect.computeBaseMass(p.getPeptide(), p.getPos(), p.getPtmmasses());
			int[] massInds = new int[base[0].length];
			for(int i= 0; i < massInds.length; i++){
				massInds[i] = (int)(prmSpect.getMassIndex(base[0][i]));
			}
	//		double totalScore = sComp.compare(t, prmSpect);
			double[] totalScore = sComp.compare(new Spectrum[]{t, t2}, new Spectrum[]{prmSpect, prmSpect2});
			double pm1 = p.getParentmass() * p.getCharge();
			double pm2 = p2.getParentmass() * p2.getCharge();
			double pm = pm1 > pm2 ? pm1 : pm2;
			MXGF mixgf = new MXGF(pm, 200, aafreq.getMassFreq());
			mixgf.setScoreHi(scores);
			mixgf.setScoreLow(scores2);
//			MXGF mixgf = new MXGF(393.1, 200, aafreq.getMassFreq());
//			mixgf.setScoreHi(getSimulatedPRM(455));
//			mixgf.setScoreLow(getSimulatedPRM2(455));
			mixgf.computePair();
//			double jProb = mixgf.getJointProb(1, 345, 393);
			double jProb = mixgf.getJointProb(totalScore[0], 
					(p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER),
					(p2.getParentmass()*p2.getCharge()-Mass.PROTON_MASS*p2.getCharge()-Mass.WATER));
			//double jProb = 0.00001;
			mixgf.computeCondPair(massInds);		
			count++;
//			//double prob = mixgf.getSpecProb(totalScore[1], (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
			double[] prob = mixgf.getSpecProb(totalScore[1], totalScore[2],
					(p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER),
					(p2.getParentmass()*p2.getCharge()-Mass.PROTON_MASS*p2.getCharge()-Mass.WATER));
			//System.out.println("score is: " + score + "\tPRM total score:\t" + totalScore);
			//double[] prob = new double[]{0,0};
			System.out.println(line + "\t" + score +"\t" + score2 + "\t" + totalScore[0] + "\t" + totalScore[1] + "\t" + totalScore[2] + "\t" + prob[0] + "\t" + prob[1] + "\t" + jProb);
			//break;
			if(count > 0){
				//break;
			}
			//return;
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		
	}
	
	public static double[] computeMXGF(Spectrum s, String peptide, String peptide2, SpectrumComparator scorer){
		return computeMXGF(s, peptide, peptide2, scorer, 0.5);
	}
	
	/**
	 * Compute spectral probability for a peptide/peptide spectrum match
	 * @param s
	 * @param peptide
	 * @param peptide2
	 * @param scorer
	 * @param tolerance
	 * @return
	 */
	public static double[] computeMXGF(Spectrum s, String peptide, String peptide2, SpectrumComparator scorer, double tolerance){
		s = new Spectrum(s);
		double resolution = 1.0;
		((SimpleProbabilisticScorer)scorer).setMinMatchedPeak(0);
		s.mergePeaks(s, 0.05);
		s.windowFilterPeaks(15, 25);
		s.computePeakRank();
		
		Peptide p = new Peptide(peptide);
		Peptide p2 = new Peptide(peptide2);
		//System.out.println("peptide is: " + p + "\t" + p.getParentmass() +"\t" + p.getCharge() + "\t" + p2 + "\t" + p2.getParentmass() + "\t" + p2.getCharge()+ "\tspectrum: " + s.parentMass + "\t" + s.charge);
		((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).combineCharge = p.getCharge() + p2.getCharge();
		((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).mode = 1;
		((SimpleProbabilisticScorer)scorer).matchTolerance  = resolution*0.5;
		
		s.parentMass = p.getParentmass();
		s.charge = p.getCharge();
		TheoreticalSpectrum t = new TheoreticalSpectrum(p);
		PRMSpectrum prmSpect = new PRMSpectrum(s, p.getCharge(), scorer, resolution);
		
		
		s.parentMass = p2.getParentmass();
		s.charge = p2.getCharge();
		TheoreticalSpectrum t1 = new TheoreticalSpectrum(p2);
		//PRMSpectrum prmSpect1 = new PRMSpectrum(s, p2.getCharge(), scorer, resolution);
		
		((SimpleProbabilisticScorer)scorer).matchTolerance  = tolerance;
		double score = scorer.compare(t, s);
		((SimpleProbabilisticScorer)scorer).matchTolerance  = resolution*0.5;
		
		((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).mode = 2;
		
		s.parentMass = p2.getParentmass();
		s.charge = p2.getCharge();
		TheoreticalSpectrum t2 = new TheoreticalSpectrum(p2);
		PRMSpectrum prmSpect2 = new PRMSpectrum(s, p2.getCharge(), scorer, resolution);
		
		
		//s = s.removeSharePeaks(t, tolerance);
		//s.computePeakRank();
		prmSpect.removeSharePeaks(s, p, p2, prmSpect.scoredSpectrum[1], prmSpect2.scoredSpectrum[1], tolerance);
		s.computePeakRank();
		prmSpect2 = new PRMSpectrum(s, p2.getCharge(), scorer, resolution);
		
		s.parentMass = p.getParentmass();
		s.charge = p.getCharge();
		TheoreticalSpectrum t3 = new TheoreticalSpectrum(p);
		//PRMSpectrum prmSpect3 = new PRMSpectrum(s, p.getCharge(), scorer, resolution);
		
		((SimpleProbabilisticScorer)scorer).matchTolerance  = tolerance;
		double score2 = scorer.compare(t2, s);
		((SimpleProbabilisticScorer)scorer).matchTolerance  = resolution*0.5;
		
		double[] scores = prmSpect.getScoredSpectrum(p.getParentmass()*p.getCharge());
		double[] scores1 = prmSpect.getScoredSpectrum(p2.getParentmass()*p2.getCharge());
		double[] scores2 = prmSpect2.getScoredSpectrum(p2.getParentmass()*p2.getCharge());
		double[] scores3 = prmSpect2.getScoredSpectrum(p.getParentmass()*p.getCharge());
		
		PRMSpectrumComparator sComp = new PRMSpectrumComparator();
		double[][] base = prmSpect.computeBaseMass(p.getPeptide(), p.getPos(), p.getPtmmasses());
		int[] massInds = new int[base[0].length];
		for(int i= 0; i < massInds.length; i++){
			massInds[i] = prmSpect.getMassIndex(base[0][i]);
		}
		
//		double totalScore = sComp.compare(t, prmSpect);
		double[] totalScore = sComp.compare(new Spectrum[]{t, t2}, new Spectrum[]{prmSpect, prmSpect2});
		double pm1 = p.getParentmass() * p.getCharge();
		double pm2 = p2.getParentmass() * p2.getCharge();
		double pm = pm1 > pm2 ? pm1 : pm2;
		MXGF mixgf = new MXGF(pm1, pm2, 200, null, resolution);
		mixgf.setScoreHi(scores);
		mixgf.setScoreLow(scores2);

		double jProb = 0.0;
		long start = (new GregorianCalendar()).getTimeInMillis();

		//mixgf.computePair();
		
		long end = (new GregorianCalendar()).getTimeInMillis();

		//jProb = mixgf.getJointProb(totalScore[0], 
		//(p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER),
		//(p2.getParentmass()*p2.getCharge()-Mass.PROTON_MASS*p2.getCharge()-Mass.WATER));
		mixgf.computeCondPair(massInds);		

	
		double[] prob = mixgf.getSpecProb(score, score2,//totalScore[1], totalScore[2],
				(p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER),
				(p2.getParentmass()*p2.getCharge()-Mass.PROTON_MASS*p2.getCharge()-Mass.WATER));
		
		mixgf = new MXGF(pm2, pm1, 200, null, resolution);
		mixgf.setScoreHi(scores1);
		mixgf.setScoreLow(scores3);
		mixgf.computeCondPair(massInds);
		
//		System.out.println("BEGIN IONS");
//		for(int i = 0; i < scores1.length; i++){
//			System.out.println(i+"\t" + (scores1[i]));
//		}
//		System.out.println("END IONS");
//		
//		System.out.println("BEGIN IONS");
//		for(int i = 0; i < scores2.length; i++){
//			System.out.println(i+"\t" + (scores2[i]));
//		}
//		System.out.println("END IONS");
		

//		double[] prob2 = mixgf.getSpecProb(totalScore[1], totalScore[2],
//				(p2.getParentmass()*p2.getCharge()-Mass.PROTON_MASS*p2.getCharge()-Mass.WATER),
//				(p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
		double[] prob2=prob;
		return new double[]{score, score2, totalScore[0], totalScore[1], totalScore[2], prob[0], prob[1], jProb, prob2[0], prob2[1], end-start};
	}
	
	protected static double[] computeMMXGF(Spectrum s, String[] peptides, SpectrumComparator scorer, double tolerance){
		s = new Spectrum(s);
		double resolution = 1.0;
		((SimpleProbabilisticScorer)scorer).setMinMatchedPeak(0);
		s.mergePeaks(s, 0.05);
		s.windowFilterPeaks(15, 25);
		s.computePeakRank();
		//s.computePeakSector();
		Peptide[] pepList = new Peptide[peptides.length];
		for(int i = 0; i < pepList.length; i++){
			pepList[i] = new Peptide(peptides[i]);
		}
		
		Peptide p = pepList[0];
		Peptide p2 = pepList[1];
		//System.out.println("peptide is: " + p + "\t" + p.getParentmass() +"\t" + p.getCharge() + "\t" + p2 + "\t" + p2.getParentmass() + "\t" + p2.getCharge()+ "\tspectrum: " + s.parentMass + "\t" + s.charge);
		((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).combineCharge = p.getCharge() + p2.getCharge();
		((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).mode = 1;
		((SimpleProbabilisticScorer)scorer).matchTolerance  = resolution*0.5;
		
		s.parentMass = p.getParentmass();
		s.charge = p.getCharge();
		TheoreticalSpectrum t = new TheoreticalSpectrum(p);
		PRMSpectrum prmSpect = new PRMSpectrum(s, p.getCharge(), scorer, resolution);

		//PRMSpectrum prmSpect1 = new PRMSpectrum(s, p2.getCharge(), scorer, resolution);
		
		((SimpleProbabilisticScorer)scorer).matchTolerance  = tolerance;
		double score = scorer.compare(t, s);
		((SimpleProbabilisticScorer)scorer).matchTolerance  = resolution*0.5;
		
		
		double[] scores = prmSpect.getScoredSpectrum(p.getParentmass()*p.getCharge());
		
		
		double pm1 = p.getParentmass() * p.getCharge();
		double pm2 = p2.getParentmass() * p2.getCharge();
		double pm = pm1 > pm2 ? pm1 : pm2;
		MXGF mixgf = new MXGF(pm1, pm2, 200, null, resolution);
		mixgf.setScoreHi(scores);
		mixgf.computeHi();
		double prob = mixgf.getSpecProb(score, (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
		
		List<Double> probsList = new ArrayList<Double>();
		probsList.add(score);
		probsList.add(prob);
		
		for(int i = 1; i < pepList.length; i++){
			Object[] prevmxgf = computeConditionalMXGF(mixgf, s, pepList[i-1], prmSpect, pepList[i], scorer, tolerance);
			prmSpect = (PRMSpectrum)prevmxgf[2];
			probsList.add((Double)prevmxgf[0]);
			probsList.add(((double[])prevmxgf[1])[1]);
		}
		double[] probs = new double[probsList.size()];
		for(int i = 0; i < probsList.size(); i++){
			probs[i] = probsList.get(i);
		}
		return probs;
		
	}
	
	protected static Object[] computeConditionalMXGF(MXGF mxgf, Spectrum s, Peptide p1, PRMSpectrum prmSpect, Peptide p2, SpectrumComparator scorer, double tolerance){
		System.out.println("Number of peaks start with: " + s.getPeak().size());
		double resolution = 1.0;
		((SimpleProbabilisticScorer)scorer).matchTolerance  = resolution/2;
		((MixturePeakScoreLearner)((SimpleProbabilisticScorer)scorer).comp).mode = 2;
		s.parentMass = p2.getParentmass();
		s.charge = p2.getCharge();
		TheoreticalSpectrum t2 = new TheoreticalSpectrum(p2);
		PRMSpectrum prmSpect2 = new PRMSpectrum(s, p2.getCharge(), scorer, resolution*2);
		
		prmSpect2.removeSharePeaks(s, p1, p2, prmSpect.scoredSpectrum[1], prmSpect2.scoredSpectrum[1], tolerance);
		s.computePeakRank();
		prmSpect2 = new PRMSpectrum(s, p2.getCharge(), scorer, resolution);
		
		
		((SimpleProbabilisticScorer)scorer).matchTolerance  = tolerance;
		double score2 = scorer.compare(t2, s);
		((SimpleProbabilisticScorer)scorer).matchTolerance  = resolution*0.5;
		
		double pm1 = p1.getParentmass() * p1.getCharge();
		double pm2 = p2.getParentmass() * p2.getCharge();
		double pm = pm1 > pm2 ? pm1 : pm2;
		MXGF mixgf = new MXGF(pm1, pm2, 200, null, resolution);
		
		double[] scores = prmSpect.getScoredSpectrum(p1.getParentmass()*p1.getCharge());
		double[] scores2 = prmSpect2.getScoredSpectrum(p2.getParentmass()*p2.getCharge());

		mixgf.setScoreHi(scores);
		mixgf.setScoreLow(scores2);
		
		double[][] base = prmSpect.computeBaseMass(p1.getPeptide(), p1.getPos(), p1.getPtmmasses());
		int[] massInds = new int[base[0].length];
		for(int i= 0; i < massInds.length; i++){
			massInds[i] = prmSpect.getMassIndex(base[0][i]);
		}
		
		mixgf.computeCondPair(massInds);		
		double[] prob = mixgf.getSpecProb(0, score2,//totalScore[1], totalScore[2],
				(p1.getParentmass()*p1.getCharge()-Mass.PROTON_MASS*p1.getCharge()-Mass.WATER),
				(p2.getParentmass()*p2.getCharge()-Mass.PROTON_MASS*p2.getCharge()-Mass.WATER));
		
		
		return new Object[]{score2, prob, prmSpect2};
	}
	
	/**
	 * Compute single-peptide spectral probability for PSM match
	 * @param s
	 * @param peptide
	 * @param scorer
	 * @return
	 */
	public static double[] computeMSGF(Spectrum s, String peptide, SpectrumComparator scorer){
		((SimpleProbabilisticScorer)scorer).setMinMatchedPeak(0);
		s.windowFilterPeaks(15, 25);
		s.computePeakRank();
		s.computePeakSector();
		if(peptide.startsWith("r")){
			peptide = peptide.substring(1);	
		}
		Peptide p = new Peptide(peptide);
		System.out.println("peptide is: " + p + "\t" + p.getParentmass() +"\t" + p.getCharge() + "\t" + "\tspectrum: " + s.parentMass + "\t" + s.charge);
		s.parentMass = p.getParentmass();
		s.charge = p.getCharge();
		TheoreticalSpectrum t = new TheoreticalSpectrum(p);
		PRMSpectrum prmSpect = new PRMSpectrum(s, p.getCharge(), scorer, 1.0);
		double score = scorer.compare(t, s);
		double[] scores = prmSpect.getScoredSpectrum(p.getParentmass()*p.getCharge());
		PRMSpectrumComparator sComp = new PRMSpectrumComparator();
		double[][] base = prmSpect.computeBaseMass(p.getPeptide(), p.getPos(), p.getPtmmasses());
		int[] massInds = new int[base[0].length];
		for(int i= 0; i < massInds.length; i++){
			massInds[i] = (int)(Math.round(base[0][i]*0.9995));
		}
		double totalScore = sComp.compare(t, prmSpect);
		double pm1 = p.getParentmass() * p.getCharge();
		MXGF mixgf = new MXGF(pm1, pm1, 200, null);
		mixgf.setScoreHi(scores);
		mixgf.setScoreLow(scores);
		mixgf.computeCondPair(massInds);		
		double prob = mixgf.getSpecProb(totalScore, (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
		System.out.println("Scan\t" + s.scanNumber + "\t" + peptide + "\t" +  + score +"\t" + totalScore + "\t" + prob);
		return new double[]{score, totalScore, prob};
	}
	
	protected static double[] computeMSGF(Spectrum s, String peptide, SpectrumComparator scorer, Spectrum scoredSpectrum ){
		if(peptide.startsWith("r")){
			peptide = peptide.substring(1);
		}
		Peptide p = new Peptide(peptide);
		System.out.println("peptide is: " + p + "\t" + p.getParentmass() +"\t" + p.getCharge() + "\t" + "\tspectrum: " + s.parentMass + "\t" + s.charge);
		s.windowFilterPeaks(6, 25);
		s.computePeakRank();
		s.parentMass = p.getParentmass();
		s.charge = p.getCharge();
		TheoreticalSpectrum t = new TheoreticalSpectrum(p);
		double pm1 = p.getParentmass() * p.getCharge();
		double[] scores = new double[(int)pm1];
		for(int i = 0; i < scoredSpectrum.getPeak().size(); i++){
			scores[i+1]=scoredSpectrum.getPeaks().get(i).getIntensity();
		}
		PRMSpectrumComparator sComp = new PRMSpectrumComparator();
		double[][] base = PRMSpectrum.computeBaseMass(p.getPeptide(), p.getPos(), p.getPtmmasses());
		int[] massInds = new int[base[0].length];
		for(int i= 0; i < massInds.length; i++){
			massInds[i] = (int)(Math.round(base[0][i]*0.9995));
		}
		double score = scorer.compare(t, s);
		double totalScore = 0;
		for(int i= 0; i < massInds.length; i++){
			if(massInds[i] > 0){
				totalScore+=scoredSpectrum.getPeak().get(massInds[i]-1).getIntensity();
			}
		}
		MXGF mixgf = new MXGF(pm1, pm1, 200, null);
		mixgf.setScoreHi(scores);
		mixgf.setScoreLow(scores);
		double jProb = 0.0;
		mixgf.computeCondPair(massInds);		
		double prob = mixgf.getSpecProb(totalScore, (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
		System.out.println("Scan\t" + s.scanNumber + "\t" + peptide + "\t" +  + score +"\t" + totalScore + "\t" + prob + "\t" + jProb);
		return new double[]{};
	}
	
	/**
	 * Test computing spectral probability for the single-peptide case
	 * This is essentially equivalent to msgf
	 * @param spectrumFile
	 * @param prmFile
	 * @param scorerFile
	 */
	public static void testComputeMSGF(String spectrumFile, String scorerFile){
		LargeSpectrumLibIterator<Spectrum> it = new LargeSpectrumLibIterator(spectrumFile);
		PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		while(it.hasNext()){
			Spectrum s = it.next();
			computeMSGF(s, s.peptide+"."+s.charge, scorer);
		}
	}
	
	protected static void testComputeMSGF2(String spectrumFile, String annotationFile, String scorerFile, int minScan, int maxScan){
		List<String> resultLines = Utils.FileIOUtils.createListFromFile(annotationFile);
		Map<Integer, String[]> table = new HashMap<Integer, String[]>();
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		int scanColumn = 1;
		int peptideColumn = 6;
		int chargeInd = 7;
		int counter = 0;
		for(int r = 0;  r < resultLines.size(); r++){
			String result = resultLines.get(r);
			if(result.startsWith("#")){
				continue;
			}
			String[] tokens = result.split("\\t");
			tokens[1] = tokens[1].replaceAll("Scan Number: ", "");
			Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[scanColumn]));
			if(s.scanNumber < minScan || s.scanNumber >= maxScan){
				continue;
			}
			//String peptide = StringUtils.getPepSeq(tokens[peptideColumn])
			//		+ "." + Integer.parseInt(tokens[chargeInd]);
			String peptide = tokens[9];
			double[] stat = computeMSGF(s, peptide, scorer);
			for(int i = 0; i < tokens.length; i++){
				System.out.print(tokens[i] + "\t");
			}
			for(int i = 0; i < stat.length; i++){
				System.out.print(stat[i] + "\t");
			}
			System.out.println();
			
			peptide = tokens[10];
			stat = computeMSGF(s, peptide, scorer);
			for(int i = 0; i < tokens.length; i++){
				System.out.print(tokens[i] + "\t");
			}
			for(int i = 0; i < stat.length; i++){
				System.out.print(stat[i] + "\t");
			}
			System.out.println();

		}

	}
	
	/**
	 * Test computing spectral probability for the single-peptide case
	 * Given the prm spectra in a file. This will produce the same probability
	 * as msgf given the prm spectrum output from it.
	 * @param spectrumFile
	 * @param prmFile
	 * @param scorerFile
	 */
	public static void testComputeMSGF(String spectrumFile, String prmFile, String scorerFile){
		LargeSpectrumLibIterator<Spectrum> it = new LargeSpectrumLibIterator(spectrumFile);
		LargeSpectrumLibIterator<Spectrum> it2 = new LargeSpectrumLibIterator(prmFile);
		PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		while(it.hasNext()){
			Spectrum s = it.next();
			Spectrum prm = it2.next();
			computeMSGF(s, s.peptide+"."+s.charge, scorer, prm);
		}
	}
	
	protected static void testSimulatedMXGF(String spectrumFile, String scorerFile, double alpha, int minCount, int maxCount){
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		PeakComparator comp = MixturePeakScoreLearner.loadComparator(scorerFile);
		//PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		List<Spectrum> specList = lib.getAllSpectrums();
		int count = 0;
		for(int i = 0; i < specList.size(); i++){
			Spectrum s1 = specList.get(i);
			for(int j = i+1; j < specList.size(); j++){
				Spectrum s2 = specList.get(j);				
				//System.out.println("count is: " + count + "\t" + minCount + "\t" + maxCount);
				if(Math.abs(s1.parentMass - s2.parentMass) < 3
						&& s1.charge == s2.charge
						&& !s1.peptide.equals(s2.peptide)){
					count++;
					if(count < minCount || count > maxCount){
						continue;
					}
					Spectrum mix = new Spectrum(s1, s2, 1.0, alpha, false);
					computeMXGF(s1, s1.peptide, s2.peptide, scorer);
					computeMXGF(s2, s1.peptide, s2.peptide, scorer);
					computeMXGF(mix, s1.peptide, s2.peptide, scorer);
				}
			}
		}
		
	}
	
	
	protected static void testRandomSimulatedMXGF(String spectrumFile, String scorerFile, double alpha, int minCount, int maxCount){
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		PeakComparator comp = MixturePeakScoreLearner.loadComparator(scorerFile);
		//PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		List<Spectrum> specList = lib.getAllSpectrums();
		int count = 0;
		while(count < maxCount){
			int i = (int)(Math.random()*specList.size()-1);
			int j = (int)(Math.random()*specList.size()-1);
			Spectrum s1 = specList.get(i);
			Spectrum s2 = specList.get(j);				
				//System.out.println("count is: " + count + "\t" + minCount + "\t" + maxCount);
				if(Math.abs(s1.parentMass - s2.parentMass) < 3
						&& s1.charge == s2.charge
						&& !s1.peptide.equals(s2.peptide)){
					count++;
					Spectrum mix10 = new Spectrum(s1, s2, 1.0, 1.0, false);
					Spectrum mix05= new Spectrum(s1, s2, 1.0, 0.5, false);
					Spectrum mix03= new Spectrum(s1, s2, 1.0, 0.3, false);
					Spectrum mix02= new Spectrum(s1, s2, 1.0, 0.2, false);
					Spectrum mix01= new Spectrum(s1, s2, 1.0, 0.1, false);
					computeMXGF(s1, s1.peptide, s2.peptide, scorer);
					computeMXGF(s2, s1.peptide, s2.peptide, scorer);
					computeMXGF(mix10, s1.peptide, s2.peptide, scorer);
					computeMXGF(mix05, s1.peptide, s2.peptide, scorer);
					computeMXGF(mix03, s1.peptide, s2.peptide, scorer);
					computeMXGF(mix02, s1.peptide, s2.peptide, scorer);
					computeMXGF(mix01, s1.peptide, s2.peptide, scorer);
				}
			}
	}
	
	/**
	 * Compute mixgf spectral probabilitiest for a set of simulated spectrum
	 * @param annotationFile
	 * @param spectrumFile
	 * @param scorerFile
	 * @param alpha
	 * @param minCount
	 * @param maxCount
	 */
	public static void testSimulatedMXGF(String annotationFile, String spectrumFile, String scorerFile, double alpha, int minCount, int maxCount){
		List<String> resultLines = Utils.FileIOUtils.createListFromFile(annotationFile);
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		PeakComparator comp = MixturePeakScoreLearner.loadComparator(scorerFile);
		//PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		List<Spectrum> specList = lib.getAllSpectrums();
		int count = 0;
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < resultLines.size(); i++){
			String line = resultLines.get(i);
			String[] tokens = line.split("\\s+");
			String peps = tokens[2] + " & " + tokens[4];
			String peptide1 = tokens[4];
			String peptide2 = tokens[2];
			String decoy1 = tokens[9];
			String decoy2 = tokens[11];
			if(lib.getSpectra(peps) == null){
				System.out.println("warning null");
				continue;
			}
			//Spectrum s1 = lib.getSpectra(tokens[2]).get(0);
			Spectrum s1 = lib.getSpectra(peps).get(0);
			count++;
			if(count < minCount || count > maxCount){
				continue;
			}
			s1.spectrumName = s1.peptide;
			double[] scores = computeMXGF(s1, peptide1, peptide2, scorer);
			double[] scores2 = computeMXGF(s1, decoy1, decoy2, scorer);
			System.out.println(line + "\t" + scores[0] +"\t" +  scores[1] + "\t" + scores[2] + "\t" + scores[3] + "\t" + scores[4] + "\t" + scores[5] + "\t" + scores[6] + "\t" + scores[7]);
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	//try to get the spectrum file from search result
	//if given a directory we automatically extract spectrumFile from search result
	//otherwise we assume the user specify the spectrum path correctly and try to use that directly
	private static String getSpectrumFromResult(String spectPath, String result){
		File f = new File(spectPath);
		if(!f.isDirectory()){
			return spectPath;
		}
		String[] tokens = result.split("\\s+");
		String spectFile = Utils.FileIOUtils.getFileName(tokens[0]);
		return spectPath+File.separator+spectFile;
	}
	/**
	 * This is the main entry point to mixgf to calculate spectral probability for mixture spectrum
	 * @param resultFile
	 * @param spectrumFile
	 * @param scorerFile
	 * @param minScan
	 * @param maxScan
	 * @param tolerance
	 * @param outfile
	 */
	public static void testComputeMXGF(String resultFile, String spectrumFile, String scorerFile, int minScan, int maxScan, double tolerance, String outfile){
		int scanIndex = 1;
		int pepIndex1 = 8;
		int pepIndex2 = 9;
		int chargeIndex1=6;
		int chargeIndex2=7;
		List<String> resultLines = Utils.FileIOUtils.createListFromFile(resultFile);
		AAFreq aafreq = null;
		PrintStream out = Utils.FileIOUtils.getOutStream(outfile);
		//nothing to do when result is empty (only header)
		if(resultLines.size() < 2){
			return;
		}
		MZXMLReader iterator = new MZXMLReader(getSpectrumFromResult(spectrumFile, resultLines.get(1)));
		
		//LargeSpectrumLibIterator specIt = new LargeSpectrumLibIterator(spectrumFile);
		//SpectrumLibMap reader = new SpectrumLibMap(spectrumFile, "MGF");
		//PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		PeakComparator comp = MixturePeakScoreLearner.loadComparator(scorerFile);
//		String trainFile = "..\\mixture_linked\\mixtures.mgf";
//		String outfile = "..\\mixture_linked\\mixtures_TOF_alpha01-10_models.o";
//		MixturePeakScoreLearner peakscorer = new MixturePeakScoreLearner(trainFile); //scorer
//		peakscorer.getMixtureIonCount();
//		PeakComparator comp = peakscorer;
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		((SimpleProbabilisticScorer)scorer).setMinMatchedPeak(0);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int count = 0;
		double time = 0;
		int specInd = 0;
		Spectrum s=null;
		for(Iterator<String> it = resultLines.iterator(); it.hasNext();){
			String line = it.next();
			String[] tokens = line.split("\\t");
			 s = iterator.getSpectrum(Integer.parseInt(tokens[scanIndex]));
			//Spectrum s = reader.getSpecByScan(Integer.parseInt(tokens[scanIndex]));
			//while(specIt.hasNext() && specInd < Integer.parseInt(tokens[scanIndex])){
			//s = (Spectrum)specIt.next();
			//	specInd++;
			//}
			
			if(s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}
			String[] peptides = new String[]{tokens[pepIndex1], tokens[pepIndex2]};
			//String[] peptides = tokens[5].split(" & ");
			//System.out.println("number of peptides " + peptides.length);
			String pep1 = tokens[pepIndex1]+"."+tokens[chargeIndex1];
			String pep2 = tokens[pepIndex2]+"."+tokens[chargeIndex2];
			//double[] scores;
			//System.out.println("peptide: " + pep1 + "\t" + pep2);
			double[] scores = computeMXGF(s, pep1, pep2, scorer, tolerance);
			time += scores[scores.length-1];
			//scores = computeMMXGF(s, peptides, scorer, tolerance);
			out.print(line +"\t");
			out.println(ArrayUtils.getString(scores, "%G"));
			count++;
			//double[] scores2 = computeMXGF(s, pep2, pep1, scorer, tolerance);
			//System.out.print(line +"\t");
			//System.out.print(Utils.StringUtils.switchMixDBStat(line) + "\t");
			//System.out.println(ArrayUtils.getString(scores2, "%G"));
			//if(count % 10 == 0){
			//	System.out.println("Computing " + count + " mix-spectral probability in " + time/1000 + " seconds");
			//}
		}
		System.out.println("Conputing  " + count+ " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static double[] getSimulatedPRM(int pm){
		double[] scores = new double[pm];
		scores[115] = 1.0;
		scores[230] = 1.0;
		scores[345] = 1.0;
		return scores;
	}
	
	public static double[] getSimulatedPRM2(int pm){
		double[] scores = new double[pm];
		scores[131] = 1.0;
		scores[262] = 1.0;
		scores[393] = 1.0;
		return scores;
	}
	
	public static void main(String[] args){
		//testComputeSingleMXGF();
		//testComputeMSGF();
//		args[0] = "../mixture_linked/MixDBv1.0/Testing/test_Kfattila_diamzxml/test.out";
//		args[1] = "../mixture_linked/MixDBv1.0/Testing/test_Kfattila_diamzxml/cyt1_400to525_b.mzXML";
//		args[2] = "../mixture_linked/yeast_simmix_alpha_generic_12_25.o";
		//args[2] = "../mixture_linked/mixtures_TOF_alpha01-10_models.o";
//		args[3] = "1";
//		args[4] = "200000";	
//		args[5] = "0.5";
		String resultFile = args[0];
		String spectrumFile = args[1];
		String scorerFile = args[2];
		int minScan = Integer.parseInt(args[3]);
		int maxScan = Integer.parseInt(args[4]);
		double tolerance = Double.parseDouble(args[5]);
		String outfile = args[6];
		testComputeMXGF(resultFile, spectrumFile, scorerFile, minScan, maxScan, tolerance, outfile);
//		args[0] = "../mixture_linked/yeast_data/klc_010908p_yeast-digest.mzXML";
//		args[0] = "../mixture_linked/yeast_simulatedmixture_samecharge_alpha1.0.mgf";
//		args[1] = "../mixture_linked/yeast_single_model_realannotated_win12_25.o";
//		args[2] = "../mixture_linked/yeast0_annotation.txt";
//		args[2] = "..\\mixture_linked\\MSGFSearches\\testprmout.txt";
//		args[1]= "../mixture_linked/yeast_simmix_alpha_generic_12_25.o";
//		args[2] = "1.0";
//		args[3] = "1"; 
//		args[4] = "10000";
//		args[5] = "../mixture_linked/testAnnotation.txt2";
		//testSimulatedMXGF(args[5], args[0], args[1], Double.parseDouble(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		//testRandomSimulatedMXGF(args[0], args[1], Double.parseDouble(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		//testComputeMSGF(args[0], args[1]);
		//testComputeMSGF2(args[1], args[0], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		//testComputeMSGF(args[0], args[2], args[1]);
	}
	
}
