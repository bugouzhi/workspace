package org.Spectrums;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
/**
 * Fully featured score learner for linked peptides.  Note there are other implementation
 * of linkedpeptide scorinig model in the package, but they are for development purpose
 * this class should be the latest implementation with more feature
 * 
 * @author Jian
 *
 */

public class LinkedPeptidePeakScoreLearner implements PeakComparator, Serializable{
	private static final long serialVersionUID = 234823048028349320L;
	private String[] ionsType = Mass.standardIonsType;
	private static int MAXRANK = 2000;
	private static int MAXLENGTH = 150;
	private static int PEPCOUNT = 2;
	private static int LINKEDTYPE = 2;
	public static final int[] DefaultRanks = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30, 35, 
			40, 45, 50, 55,60,65, 70,75,80, 85, 90,100,110, 120, 130, 140, 150,MAXRANK};
	private int[] rankInterval = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30, 35, 
			40, 45, 50, 55,60,65, 70,75,80, 85, 90,100,110, 120, 130, 140, 150,MAXRANK};//160,170,180,190,200,210,220,230,240,250, MAXRANK};
	private int[] lengthInterval = {1,MAXLENGTH};
	//to define an convention, here mass error means theoretical mass - actual mass of a peak
	private double[] massToleranceInterval = {-0.51,0.51}; //need to setup this automatically in the future
	private double[] massErrorInterval = {-0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55};
	//private double[] massErrorInterval = {-0.055, -0.045, -0.035, -0.025, -0.015, -0.005, 0.005, 0.015, 0.025, 0.035, 0.045, 0.055};
	private HashMap<String, Integer> ionIndex;
	private int commonMinCharge = 1;
	private int commonMaxCharge = 6;
	private int linkedMinCharge = 2;
	private int linkedMaxCharge = 8;
	private Map<String, Double> priorProbs;
	private LookUpTable table;
	private LookUpTable errorModel;
	private SpectrumLib annotatedSet;
	private String annnotateSetFile;
	public static boolean offSetFreq = false; //toggle whether we are learning the ion stat in discovery mode
	public boolean generalLinkedMode=false;  //to simulate general linked model rather than the SUMO specific model, where we turn off the SUMO tag model, only use the peptide scoring model
	public int peptideMode = 0;
	
	public LinkedPeptidePeakScoreLearner(String file){
		this.annnotateSetFile = file;
		initialize();
	}
	
	public LinkedPeptidePeakScoreLearner(SpectrumLib lib){
		this.annotatedSet = lib;
		initialize();
	}
	
	
	private void initialize(){
		this.initializeIonIndexTable();
		this.table = initializeTable();
		this.initializeErrorModel();
	}
	
	/**
	 * Initialize the fragment mass tolerance error models
	 */
	private void initializeErrorModel(){
		this.errorModel = new LookUpTable(
			new int[] {	this.rankInterval.length, this.massErrorInterval.length});    // assumed error model indenpendent of ranks
	}
	
	/**
	 * Initialize the fragment ion type models
	 */
	private void initializeIonIndexTable(){
		this.ionIndex = new HashMap<String, Integer>();
	//	this.ionsType = new String[Mass.standardIonsType.length+1];
		for(int i = 0; i < ionsType.length; i++){
			ionIndex.put(ionsType[i], new Integer(i));
			//System.out.println("storing ion: " + ionsType[i]);
		}
	}
		
	private LookUpTable initializeTable(){
		LookUpTable table = new LookUpTable( 
			new int[] {this.PEPCOUNT, this.LINKEDTYPE, this.linkedMaxCharge, lengthInterval.length, this.commonMaxCharge,  
					this.ionsType.length+1, rankInterval.length+1, this.massToleranceInterval.length}); //extra slot at the beginning for noise model
		return table;
	}
	
	@Override
	/**
	 * Score when comparing a theoretical peak to an observed peak
	 */
	public double compare(Peak p1, Peak p2) {
		if(p1 instanceof LabelledPeak && !(p1 instanceof MixturePeak)){
			return compareSingle(p1,p2);
		}
		MixturePeak lp = (MixturePeak)p1;
		//System.out.println("comparing score: " + p1 + "\t" + p2);
		if(p1 == null ){
			return 0;
		}else{
			int[] index =  getIndex(lp, p2);
			int[] index2 = getNoiseIndex(lp.getPep(), p2, lp.getParent().charge, lp.getPeptideIndex()); //assuem if we were to match the peak to noise
			int[] errorIndex = getErrorIndex(lp, p2);
//			int[] index2 = getNoiseIndex(lp, p2); //assuem if we were to match the peak to noise
			double score = this.table.get(index);
			double score2 = this.table.get(index2);
			if(p2 != null){
				double errorScore = this.errorModel.get(errorIndex);
				//System.out.print("peptide " + lp.getPeptideIndex()  +  "\t" +  lp +  "\t~\t" + p2 + "\tscore: " + score + "\t" + score2 + "\t" + errorScore);
				score *= errorScore;
				if(score == 0){
					score = 0.00001;
					//return 0;
				}
				
				score2 = score2 / (this.massErrorInterval.length-1);
				//System.out.println("\t" +  (1.0/(this.massErrorInterval.length-1)) + "\t" + Math.log(score/score2));
			}
//			System.out.println("score: " + score);
//			System.out.println("score2: " + score2);
			if(Double.isNaN(score)){
				return 0;
				//score = 0.00001;
			}
			if(score2 == 0){
				//System.out.println("warning noise score is zero");
				//System.out.println("index is: " + lp.getParent().charge + "\t" + lp.getPeptideIndex() + "\t"
				//		+ lp.getPep() + "\t" + lp.getType() + "\t" + lp.getCharge());
				score2 = 0.00001;
			}
			if(score < score2){
				//System.out.println("score is: " + score + " score2 is: " + score2);
				//return 0.0;
			}

			//we need to take care of cases, if peaks is very unlikely, consider it not matched rather than match some 
			//very bad peaks
			double ratio = Math.log(score/score2);
			int[] indexMissing = getIndex(lp, null);
			double scoreMissing = this.table.get(indexMissing);
			int[] indexMissing2 = getNoiseIndex(lp.getPep(), null, lp.getParent().charge, lp.getPeptideIndex());
			double scoreMissing2 = this.table.get(indexMissing2);
			double ratio2 = Math.log(scoreMissing/scoreMissing2);
			//System.out.println("\t" + ratio +"\t" + ratio2);
			if(ratio2 > ratio){
				return ratio2; 
			}else{
				return ratio;
			}
		}
	}
	
	/**
	 * This handle the case when we have only "half" of the linked peptide at hand.  This will
	 * score the peak base on the "first-peptide" model (i.e. in the mixture framework it is
	 * the higher abundance peptide and in the SUMO framework this is the model for the peptide
	 * i.e. not the SUMO tag model)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public double compareSingle(Peak p1, Peak p2) {
		LabelledPeak lp = (LabelledPeak)p1;
		//System.out.println("comparing score: " + p1 + "\t" + p2);
		if(p1 == null ){
			return 0;
		}else{
			int[] index =  getIndex(lp, p2);
			int[] index2 = getNoiseIndex(lp.getPep(), p2, lp.getPep().getCharge(), 0); //assuem if we were to match the peak to noise
			int[] errorIndex = getErrorIndex(lp, p2);
//			int[] index2 = getNoiseIndex(lp, p2); //assuem if we were to match the peak to noise
			double score = this.table.get(index);
			double score2 = this.table.get(index2);
			if(p2 != null){
				//System.out.println("peptide " + lp +  "\trank: " + p2.getRank() + " score: " + score + "\t" + score2 + "\t" + Math.log(score/score2));
				double errorScore = this.errorModel.get(errorIndex);
				score *= errorScore;
				if(score == 0){
					score = 0.00001;
				}
				score2 = score2 / (this.massErrorInterval.length-1);
				
				//System.out.println("peptide " + 0  +  "\t" +  lp +  "\trank: " + p2.getRank() + " score: " + score + "\t" + score2 + "\t" + Math.log(score/score2) + "\t" + errorScore);
			}
//			System.out.println("score: " + score);
//			System.out.println("score2: " + score2);
			if(Double.isNaN(score)){
				return 0;
				//score = 0.00001;
			}
			if(score2 == 0){
				//System.out.println("warning noise score is zero");
				score2 = 0.00001;
			}
			if(score < score2){
				//System.out.println("score is: " + score + " score2 is: " + score2);
				//return 0.0;
			}
//			if(p2 != null){
//				System.out.println("scoring: " + lp + "\t" + p2 + "\t:\t" + Math.log(score/score2));
//			}else{
//				System.out.println("scoring: " + lp +  ":\t" + Math.log(score/score2));
//			}
			return Math.log(score/score2);
		}
	}
	
	
	/**
	 * Given a theoretical fragment ions, get the index for the ion type
	 * @param lp
	 * @return
	 */
	protected int getIonIndex(LabelledPeak lp){
		if(!this.ionIndex.containsKey(lp.getType())){
			throw new IllegalArgumentException("Invalide ion type " + lp.getType());
		}
		return ionIndex.get(lp.getType()).intValue();
	}
	
	/**
	 * Peptide separated into different type based on length of the peptide
	 * @param p
	 * @return
	 */
	protected int getPeptideLength(Peptide p){
			return ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
	}
	
	/**
	 * return the value in the high-dimensional table given a index
	 * @param index
	 * @return
	 */
	protected double getValue(int[] index){
		return this.table.get(index);
	}
	//scoring table dimesnion
	//peptide index
	//linked type (i.e. whether it is linked peaks or not)
	//peptide charge
	//length
	//peakcharge
	//ion type
	//mass error
	
	/**
	 * Given a theoretical peak and a real peak pair, return its corresponding index in the score table
	 * @param lp
	 * @param realPeak
	 * @return
	 */
	private int[] getIndex(MixturePeak lp, Peak realPeak){
		int rankIndex;
		int errorIndex;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
			errorIndex = ArrayUtils.getIntervalIndex(
					lp.getMass() - realPeak.getMass(), this.massToleranceInterval); 
		}
		if(lp == null){
			return new int[]{0,0,0,0,rankIndex};
		}
		Peptide p = lp.getPep();
		int peptideLength = getPeptideLength(lp.getPep());//ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		///System.out.println("peptide is: " + p.getPeptide() + " index : " + peptideLength);
		int isLinked = 0;
		if(TheoreticalSpectrum.isLinkedPeak(p, lp)){
			isLinked = 1;
			//System.out.println("peak is linked " + lp + "linked position " + lp.getPep().getLinkedPos());
		}
		int peptideCharge = lp.getPep().getCharge()-1;
		if(lp.getPep().getCharge() > 5){
			peptideCharge=4;
		}
//		if(lp.getPep().getCharge() < 2){
//			peptideCharge=2;
//		}
		int peakCharge = lp.getCharge()-1;
		if(lp.getCharge() > 5){
			peakCharge = 4;
		}
		int ionIndex = getIonIndex(lp)+1;
		//System.out.println("peak " + lp.getPeptideIndex() +"\t" + isLinked + "\t" + peptideCharge +"\t" + peptideLength +"\t" + peakCharge +"\t" + ionIndex +"\t" + rankIndex +"\t" + errorIndex);
		if(this.generalLinkedMode){
			return new int[]{0, isLinked, peptideCharge, peptideLength, peakCharge, ionIndex, rankIndex, errorIndex};
		}
		return new int[]{lp.getPeptideIndex(), isLinked, peptideCharge, peptideLength, peakCharge, ionIndex, rankIndex, errorIndex};
	}
	
	/**
	 * Same as above but do for the case when there is only half of the linked peptide (so LabelledPeak is input instead of MixturePeak)
	 * @param lp
	 * @param realPeak
	 * @return
	 */
	private int[] getIndex(LabelledPeak lp, Peak realPeak){
		int rankIndex;
		int errorIndex;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
			errorIndex = ArrayUtils.getIntervalIndex(
					lp.getMass() - realPeak.getMass(), this.massToleranceInterval); 
		}
		if(lp == null){
			return new int[]{0,0,0,0,rankIndex};
		}
		Peptide p = lp.getPep();
		int peptideLength = getPeptideLength(lp.getPep());//ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		///System.out.println("peptide is: " + p.getPeptide() + " index : " + peptideLength);
		int isLinked = 0;
		//System.out.println("checking peak " + lp);
		if(TheoreticalSpectrum.isLinkedPeak(p, lp)){
			isLinked = 1;
			//System.out.println("peak is linked " + lp + "linked position " + lp.getPep().getLinkedPos());
		}
		
		int peptideCharge = lp.getPep().getCharge()-1;
		if(lp.getPep().getCharge() > 5){
			peptideCharge=4;
		}
//		if(lp.getPep().getCharge() < 2){
//			peptideCharge=2;
//		}
		int peakCharge = lp.getCharge()-1;
		if(lp.getCharge() > 5){
			peakCharge = 4;
		}
		int ionIndex = getIonIndex(lp)+1;
		//System.out.println("peak " + 0 +"\t" + isLinked + "\t" + peptideCharge +"\t" + peptideLength +"\t" + peakCharge +"\t" + ionIndex +"\t" + rankIndex +"\t" + errorIndex);
		return new int[]{this.peptideMode, isLinked, peptideCharge, peptideLength, peakCharge, ionIndex, rankIndex, errorIndex};
	}
	
	/**
	 * Get the index in the scoring table given the error between matched theoretical and real peak
	 * @param lp
	 * @param realPeak
	 * @return
	 */
	private int[] getErrorIndex(LabelledPeak lp, Peak realPeak){
		int rankIndex;
		int errorIndex;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
			errorIndex = ArrayUtils.getIntervalIndex(
					lp.getMass() - realPeak.getMass(), this.massErrorInterval); 

		}
		return new int[]{rankIndex, errorIndex};
	}
	
	/**
	 * Given a  observed peak, return the noise score, which is the null
	 * model that the observed peak is generated from noise
	 * @param p
	 * @param realPeak
	 * @param combineCharge
	 * @param peptideIndex
	 * @return
	 */
	private int[] getNoiseIndex(Peptide p, Peak realPeak, int combineCharge, int peptideIndex){
		int rankIndex;
		int errorIndex;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
		}
		int peptideLength = getPeptideLength(p);//ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		int peptideCharge = p.getCharge()-1;
		if(p.getCharge() > 4){
			peptideCharge = 3;
		}
		int peakCharge = 0;
		int ionIndex = 0;
		return new int[]{peptideIndex, 0, peptideCharge, peptideLength, peakCharge, ionIndex, rankIndex, 0};
	}
	
	/**
	 * Noise model for mixture model
	 * @param lp
	 * @param realPeak
	 * @param combineCharge
	 * @param peptideIndex
	 * @return
	 */
	private int[] getNoiseIndex(MixturePeak lp, Peak realPeak, int combineCharge, int peptideIndex){
		int rankIndex;
		int errorIndex;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
			errorIndex = ArrayUtils.getIntervalIndex(
					lp.getMass() - realPeak.getMass(), this.massToleranceInterval); 

		}
		if(lp == null){
			return new int[]{0,0,0,0,0,rankIndex};
		}
		Peptide p = lp.getPep();
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		///System.out.println("peptide is: " + p.getPeptide() + " index : " + peptideLength);
		int peptideCharge = lp.getPep().getCharge()-1;
		int peakCharge = lp.getCharge()-1;
		int ionIndex = 0;
		return new int[]{lp.getPeptideIndex(), lp.getParent().charge, peptideLength, peptideCharge, peakCharge, ionIndex, rankIndex, errorIndex};
	}
	
	private Iterator<Spectrum> getAnnotatedIterator(){
		if(this.annotatedSet != null){
			return this.annotatedSet.getSpectrumList().iterator();
		}else{
			return new LargeSpectrumLibIterator(this.annnotateSetFile);
		}
	}
	
	//This is the entry method to the learner, it learn the statistics from
	//the set of annotated spectra used to initialize the model
	public void getLinkedIonCount(){
		LookUpTable totalCount = initializeTable();
		int count = 0;
		for(Iterator it = getAnnotatedIterator(); it.hasNext();){
			Spectrum s = (Spectrum)it.next();
//			if(s.modMass > 0){
//				continue;
//			}
			//Mass.DSSLINKER_MASS = -116.0430;
			Mass.DSSLINKER_MASS = -1*Mass.WATER;
			s.removePrecursors(0.5);
			s.windowFilterPeaks(10, 25);
			//s = new DeconvolutedSpectrum(s);
			//TheoreticalSpectrum.deconvolutedMode=true;
			s.computePeakRank();
			//s.peptide = s.peptide+"--GG";        //ubiq
			
			s.peptide = s.peptide.replace(" & ", "--");
			String[] peps = s.peptide.split("--");
			
			//System.out.println("peptide is: " + s.peptide);
			
			//need to creat linkedpeptides
			if(!s.peptide.contains("7") || s.charge > 5){
				continue;
			}
			int position1 = 1;
			int position2 = 6;
			if(s.peptide.contains("LRAK")){
				position1 = 1;
				System.out.println("lib_sumo1");
			}else if(s.peptide.contains("TALH")){
				position1 = 6;
				System.out.println("lib_sumo2");
			}else if(s.peptide.contains("FRAK")){
				position1 = 3;
				System.out.println("lib_sumo3");
			}else{
				position1 = -1;
			}
			
			Peptide p = new Peptide(peps[0],1);                                 //ubiq
//			double[] ptmmasses = p.getPtmmasses();
//			for(int i = 0 ; i < ptmmasses.length; i++){
//				if(ptmmasses[i] == 114.043 || ptmmasses[i] == 122.084){
//					position1 = p.getPos()[i];
//					System.out.println("peptide: " + peps[0] + " position: " + position1);
//				}
//			}
//			position2=2;                       //ubiq
//			
//			
//			for(int i = 0; i < p.getPeptide().length(); i++){
//				if(p.getPeptide().charAt(i) == 'C'){
//					position1=i+1;
//					System.out.println("peptide: " + peps[0] + " position: " + position1);
//					break;
//				}
//			}
//			
//			Peptide p2 = new Peptide(peps[1],1);                                 //disulfide
//			for(int i = 0; i < p2.getPeptide().length(); i++){
//				if(p2.getPeptide().charAt(i) == 'C'){
//					position2=i+1;
//					System.out.println("peptide: " + peps[1] + " position: " + position2);
//					break;
//				}
//			}
			
			if(s.charge > 6 || position1 < 0){
				continue;
			}
			System.out.println("peptide is: " + s.peptide);
			LinkedPeptide linked = new LinkedPeptide(s.peptide, s.charge, position1, position2);
			System.out.println("peptide is: " + linked + "\t" + s.peptide);
			System.out.println("peptide1 is: " + linked.peptides[0]);
			System.out.println("peptide2 is: " + linked.peptides[1]);
			//LinkedPeptide linked = new LinkedPeptide(peptides[0]+"--Z", s.charge, position1, 1);
			TheoreticalSpectrum t = new TheoreticalSpectrum(linked.peptides[0], linked.peptides[1], (short)s.charge, false, Mass.DSSDANGLE_MASS);
			SimpleMatchingGraph matchingG = t.getMatchGraph(s, 0.5);
			//removeSUMO(matchingG, s);
			//matchingG = t.getMatchGraph(s, 0.05);
			this.getIonsCount(matchingG, peps[0], peps[1]);
			System.out.println(s.scanNumber + "\t" + s.parentMass + s.charge);
			String[] peptides = s.peptide.split("--");
			//t.analyzeMixtureAnnotation(s, peptides[0],  peptides[1]);	
			//we counted each peptide as first and second each time for a more generalized models
			//TheoreticalSpectrum t2 = new TheoreticalSpectrum(linked.peptides[1], linked.peptides[0], (short)s.charge, false);
			//SimpleMatchingGraph matchingG2 = t2.getMatchGraph(s, 0.5);
			//this.getIonsCount(matchingG2, peps[1], peps[0]);
			
			count++;
			if(count % 1000 == 0){
				System.out.println("Finish Analyzing " + count);
			}
			if(count == 20000){
				//break;
			}
			//return;
		}
		this.normalizeCount();
		this.normalizeErrorModel();
		printIonTable();
		printErrorTable();
	}
	
	public void removeSUMO(SimpleMatchingGraph matchingG, Spectrum s){
		List<Peak> toBeRemoved = new ArrayList();
		for(Iterator it = matchingG.vertexSet(SimpleMatchingGraph.Observed).iterator(); 
			it.hasNext();){
			Peak p = (Peak)it.next();
			int type = 0;
			for(Iterator it2 = matchingG.getNeighbors(p).iterator(); it2.hasNext();){
				MixturePeak neigh = (MixturePeak)it2.next();
				//System.out.println("index: " + neigh.getPeptideIndex());
				if(neigh.getPeptideIndex() == 0){
					if(type == 0){
						type=1;
					}
					if(type==2){
						type=3;
					}
				}
				if(neigh.getPeptideIndex() == 1){
					if(type == 0){
						type=2;
					}
					if(type==1){
						type=3;
					}
				}
			}
			if(type == 2){
				toBeRemoved.add(p);
			}
		}
		System.out.println("before: " + s.getPeak().size());
		s.getPeak().removeAll(toBeRemoved);
		System.out.println("after: " + s.getPeak().size());
		s.computePeakRank();
	}
	
	/**
	 * Normalized the frequency count so each entry in the score table is a proper probability
	 */
	private void normalizeCount(){
		for(int peptide = 0; peptide < this.PEPCOUNT; peptide++){
			for(int linkedType = 0; linkedType < this.LINKEDTYPE; linkedType++){
				for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
					for(int length = 0; length < this.lengthInterval.length; length++){
						for(int peakCharge = 0; peakCharge < this.commonMaxCharge; peakCharge++){
							for(int ionIndex = 0; ionIndex < this.ionsType.length+1; ionIndex++){
								double sum = 0.0;
								for(int rank = 0; rank < this.rankInterval.length+1; rank++){
									for(int noise = 0; noise < this.massToleranceInterval.length; noise++){
										int[] index = {peptide, linkedType, pepCharge, length, peakCharge, ionIndex, rank, noise};
										sum += this.table.get(index);
									}
								}
								for(int rank = 0; rank < this.rankInterval.length+1; rank++){
									for(int noise = 0; noise < this.massToleranceInterval.length; noise++){
										int[] index = {peptide, linkedType, pepCharge, length, peakCharge, ionIndex, rank, noise};
										double count = this.table.get(index);
										if(count > sum) System.out.println("waring sum is smaller: " + count + " out of " + sum);
										this.table.put(index, count/sum);
									}
								}

							}
						}
					}
				}
			}
		}
		
		//normalize for noises
		for(int peptide = 0; peptide < this.PEPCOUNT; peptide++){
			for(int linkedType = 0; linkedType < this.LINKEDTYPE; linkedType++){
				for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
					for(int length = 0; length < this.lengthInterval.length; length++){
						double sum = 0.0;
						for(int rank = 0; rank < this.rankInterval.length+1; rank++){
							int[] index = {peptide, linkedType, pepCharge, length, 0, 0, rank, 0};
							sum += this.table.get(index);
						}
						sum /= (1-0.93);
						int[] noiseInd = {peptide, linkedType, pepCharge, length, 0, 0, 0,0};
						this.table.put(noiseInd, sum*0.95);
						for(int rank = 0; rank < this.rankInterval.length+1; rank++){
							int[] index = {peptide, linkedType, pepCharge, length, 0, 0, rank, 0};
							double count = this.table.get(index);
							this.table.put(index, count/sum);
						}
					}
				}
			}
		}
	}
	
	/**
	 * Get the ion statistics for a mixture spectrum match
	 * @param g
	 * @param pep1
	 * @param pep2
	 */
	private void getIonsCount(SimpleMatchingGraph g, String pep1, String pep2){
		Set vertices = g.vertexSet(2);
		Iterator it = vertices.iterator();
		Peak p, realPeak;
		MixturePeak lp;
		Integer c, one = new Integer(1);
		Peptide pep = null;
		//String[] peps1 = pep1.split("\\.");
		//String[] peps2 = pep2.split("\\.");
		//int combineCharge = Integer.parseInt(peps1[1]) 
		//	+ Integer.parseInt(peps2[1]);
		int peptideIndex = 0;
		
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (MixturePeak)p;
				pep = lp.getPep();
				Set<Peak> neighbors = g.getNeighborSet(lp);
				if(lp.getPep().getPeptide().equals(pep1)){
					peptideIndex = 0;
				}else{
					peptideIndex = 1;
				}
				if(neighbors.size() == 0){
					int[] index = this.getIndex(lp, null);
					this.table.incrementIonCount(index);
				}
			}
		}
		
		it = g.vertexSet(1).iterator();
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p.getRank() > 100 && offSetFreq){
				continue;
			}
			Set<Peak> neighbors = g.getNeighborSet(p);
			if(neighbors.size() == 0){				
				//int[] index = this.getIndex(null, p);
				int[] index = this.getNoiseIndex(pep, p, 0, 0);
				int[] index2 = this.getNoiseIndex(pep, p, 0, 1);

				this.table.incrementIonCount(index); //we need to count noise twice one for each peptide
				this.table.incrementIonCount(index2);
			}else{
				double massDiff=0, min=1000;
				Peak closest = null;
				for(Iterator<Peak> iter = neighbors.iterator(); iter.hasNext();){
					Peak neigh = iter.next();
					//int[] index = this.getIndex((LabelledPeak)neigh, p);
					//this.table.incrementIonCount(index, 1/neighbors.size());
					massDiff = Math.abs(neigh.getMass() - p.getMass());
					closest = massDiff < min ? neigh : closest;
					min = massDiff < min ? massDiff : min;
					if(offSetFreq){
						if(((LabelledPeak)neigh).getPep().getPeptide().equals(pep1)){
							peptideIndex = 0;
						}else{
							peptideIndex = 1;
						}
						int[] index = this.getIndex((MixturePeak)neigh, p);
						this.table.incrementIonCount(index);
						int[] errorIndex = this.getErrorIndex((LabelledPeak)neigh, p);
						this.errorModel.incrementIonCount(errorIndex);
					}
				}
				if(!offSetFreq){
					if(((LabelledPeak)closest).getPep().getPeptide().equals(pep1)){
						peptideIndex = 0;	
					}else{
						peptideIndex = 1;
					}
					int[] index = this.getIndex((MixturePeak)closest, p);
					this.table.incrementIonCount(index);
					int[] errorIndex = this.getErrorIndex((LabelledPeak)closest, p);
					this.errorModel.incrementIonCount(errorIndex);
				}
			}
		}
		//building noise model		
	}
	
	/**
	 * Normalize the statistics in the errror model
	 */
	private void normalizeErrorModel(){
			for(int rankIndex = 0; rankIndex < this.rankInterval.length; rankIndex++){
				double sum = 0.0;
				for(int error = 0; error < this.massErrorInterval.length; error++){
					int[] index = {rankIndex, error};
						sum += this.errorModel.get(index);
				}
				for(int error = 0; error < this.massErrorInterval.length; error++){
					int[] index = {rankIndex, error};
					double count = this.errorModel.get(index);
					if(count > sum) System.out.println("waring sum is smaller: " + count + " out of " + sum);
					this.errorModel.put(index, count/sum);
				}
			}
	}
	
	private void printIonTable(){
		for(int peptide = 0; peptide < this.PEPCOUNT; peptide++){
			for(int linkedType = 0; linkedType < this.LINKEDTYPE; linkedType++){
				for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
					for(int length = 0; length < this.lengthInterval.length-1; length++){
						for(int peakCharge = 0; peakCharge < this.commonMaxCharge; peakCharge++){
							for(int ionIndex = 0; ionIndex < this.ionsType.length; ionIndex++){
								double total = 0.0;
								String label;
								if(peptide == 0){
									label = "Peptide1";
								}else{
									label = "Peptide2";
								}
								
								if(linkedType == 0){
									label = label + "(regular)";
								}else{
									label = label + "(linked)";
								}
								StringBuffer header = new StringBuffer(label+ this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1));
								StringBuffer prob = new StringBuffer(header + " prob:\t");
								StringBuffer background = new StringBuffer(header + " noise:\t");
								for(int rank = 0; rank < this.rankInterval.length; rank++){
									for(int noise = 0; noise < this.massToleranceInterval.length-1; noise++){
										int[] index = {peptide, linkedType, pepCharge, length, peakCharge, ionIndex+1, rank+1, noise};
										prob.append(this.table.get(index) + "\t");
										background.append(this.table.get(new int[]{peptide, 0, pepCharge, length, 0, 0, rank+1, 0}) + "\t");
										total += this.table.get(index);
									}
								}
								System.out.println(prob.toString());
								System.out.println(background.toString());
								System.out.println(label+ this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1) + 
										" total:  "      + total); 
							}
						}
					}
				}
			}
		}
	}
	
	private void printErrorTable(){
		for(int rankIndex = 0; rankIndex < this.rankInterval.length; rankIndex++){
			for(int error = 0; error < this.massErrorInterval.length; error++){
				int[] index = {rankIndex, error};
				System.out.println("rank" + rankIndex +  
				" error: " + error + ": "    
				+ this.errorModel.get(index));
			}
		}
	
	}
	
	public void writeLibToFile(String outfile){
		try{
			BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(outfile));
			ObjectOutputStream oo = new ObjectOutputStream(bo);
		    oo.writeObject(this);
		    oo.flush();
		    oo.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		}
	}
	
	/**
	 * Load the scoring model from a file (that is previous created)
	 * @param file
	 * @return
	 */
	public static LinkedPeptidePeakScoreLearner loadComparator(String file){
		try{
			BufferedInputStream bi = new BufferedInputStream(new FileInputStream(file));
			ObjectInputStream oi = new ObjectInputStream(bi);
		    Object o = oi.readObject();
		    return (LinkedPeptidePeakScoreLearner)o;
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	public int[] getRankInterval() {
		return this.rankInterval;
	}
	
	/**
	 * Get the score for a pair of matched theoretical and observed peak
	 * same as the compare method but not using the PeakComparator interface
	 * @param p1
	 * @param p2
	 * @return
	 */
	public double getScore(Peak p1, Peak p2) {
		if(p1 instanceof LabelledPeak && !(p1 instanceof MixturePeak)){
			return compareSingle(p1,p2);
		}
		MixturePeak lp = (MixturePeak)p1;
		//System.out.println("comparing score: " + p1 + "\t" + p2);
		if(p1 == null ){
			return 0;
		}else{
			int[] index =  getIndex(lp, p2);
			int[] index2 = getNoiseIndex(lp.getPep(), p2, lp.getParent().charge, lp.getPeptideIndex()); //assuem if we were to match the peak to noise
//			int[] errorIndex = getErrorIndex(lp, p2);
//			int[] index2 = getNoiseIndex(lp, p2); //assuem if we were to match the peak to noise
			double score = this.table.get(index);
			double score2 = this.table.get(index2);
			if(p2 != null){
//				double errorScore = this.errorModel.get(errorIndex);
				//System.out.print("peptide " + lp.getPeptideIndex()  +  "\t" +  lp +  "\t" + p2 + "\tscore: " + score + "\t" + score2);
//				score *= errorScore;
				if(score == 0){
					score = 0.00001;
					//return 0;
				}
//				score2 = score2 / (this.massErrorInterval.length-1);
				//System.out.println("\t" +   + Math.log(score/score2));
			}
//			System.out.println("score: " + score);
//			System.out.println("score2: " + score2);
			if(Double.isNaN(score)){
				return 0;
				//score = 0.00001;
			}
			if(score2 == 0){
				//System.out.println("warning noise score is zero");
				//System.out.println("index is: " + lp.getParent().charge + "\t" + lp.getPeptideIndex() + "\t"
				//		+ lp.getPep() + "\t" + lp.getType() + "\t" + lp.getCharge());
				score2 = 0.00001;
			}
			if(score < score2){
				//System.out.println("score is: " + score + " score2 is: " + score2);
				//return 0.0;
			}

			//we need to take care of cases, if peaks is very unlikely, consider it not matched rather than match some 
			//very bad peaks
			double ratio = Math.log(score/score2);
			int[] indexMissing = getIndex(lp, null);
			double scoreMissing = this.table.get(indexMissing);
			int[] indexMissing2 = getNoiseIndex(lp.getPep(), null, lp.getParent().charge, lp.getPeptideIndex());
			double scoreMissing2 = this.table.get(indexMissing2);
			double ratio2 = Math.log(scoreMissing/scoreMissing2);
			if(ratio2 > ratio){
				return ratio2; 
			}else{
				return ratio;
			}
			//return ratio;
		}
	}

	/**
	 * Get the bins for the mass error model
	 * @return
	 */
	public double[] getMassErrorInterval() {
		return this.massErrorInterval;
	}

	/**
	 * Get score from the for fragment mass error model
	 * @param p1
	 * @param p2
	 * @return
	 */
	public double getErrorScore(Peak p1, Peak p2) {
		if(p1 instanceof LabelledPeak && !(p1 instanceof MixturePeak)){
			return compareSingle(p1,p2);
		}
		MixturePeak lp = (MixturePeak)p1;
		//System.out.println("comparing score: " + p1 + "\t" + p2);
		if(p1 == null ){
			return 0;
		}else{
			int[] errorIndex = getErrorIndex(lp, p2);
			if(p2 != null){
				double errorScore = this.errorModel.get(errorIndex);
				//System.out.print("peptide " + lp.getPeptideIndex()  +  "\t" +  p2 + "\t" + errorScore);
				double errorScore2 = 1.0 / (this.massErrorInterval.length-1);
				//System.out.println("\t" +  (1.0/(this.massErrorInterval.length-1)) + "\t" + Math.log(errorScore/errorScore2));
				return Math.log(errorScore/errorScore2);
			}
		}
		return 0.0;
	}
	
	/**
	 * This method create a linked-peptide scoring model from annotated MGF and 
	 * print the scoring model in a file
	 */
	public static void getLinkedModel(){
		TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		MixtureSpectrumScorer scorer2 = (MixtureSpectrumScorer)SpectrumUtil.getLinkedPeptideScorer("..\\mixture_linked\\lib_sumo_spectra_training_part1.mgf");
		((LinkedPeptidePeakScoreLearner)scorer2.comp).writeLibToFile("..\\mixture_linked\\test_model.o");
	}
	
}
