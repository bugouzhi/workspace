package org.Spectrums;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.graph.SimpleGraph;

import Utils.FileIOUtils;

/**
 * Similar to SpectrumIonRankLearner but reimplemented using
 * Lookup table
 * @author jian wang
 *
 */
public class RankBaseScoreLearner implements PeakComparator, Serializable{
	private static final long serialVersionUID = 7432728402384023L;
	private String[] ionsType = Mass.standardIonsType;
	private static int MAXRANK = 2000;
	private static int MAXLENGTH = 150;
	//private int[] rankInterval = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,60,100,MAXRANK};
	private int[] rankInterval = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30, 35, 
		40, 45, 50, 55,60,65, 70,75,80, 85, 90,100,110, 120, 130, 140, 150,	MAXRANK};
	private int[] lengthInterval = {1,12,MAXLENGTH};
	//to define an convention, here mass error means theoretical mass - actual mass of a peak
	private double[] massToleranceInterval = {-0.51,0.51}; //need to setup this automatically in the future
	private double[] massErrorInterval = {-0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55};
	//private double[] massErrorInterval = {-0.55, -0.045, -0.035, -0.025, -0.015, -0.005, 0.005, 0.015, 0.025, 0.035, 0.045, 0.55};
	private double[] sector = {0,1,2};
	public double[] getMassErrorInterval() {
		return massErrorInterval;
	}

	public void setMassErrorInterval(double[] massErrorInterval) {
		this.massErrorInterval = massErrorInterval;
	}

	private HashMap<String, Integer> ionIndex;
	private int commonMinCharge = 1;
	private int commonMaxCharge = 5;
	private int linkedMinCharge = 2;
	private int linkedMaxCharge = 6;
	private Map<String, Double> priorProbs;
	private LookUpTable table;
	private LookUpTable errorModel;
	public LookUpTable getErrorModel() {
		return errorModel;
	}

	public void setErrorModel(LookUpTable errorModel) {
		this.errorModel = errorModel;
	}

	private SpectrumLib annotatedSet;
	private String annnotateSetFile;
	public RankBaseScoreLearner(SpectrumLib lib){
		this.annotatedSet = lib;
		initialize();
	}
	
	public RankBaseScoreLearner(String specLibFile){
		this.annnotateSetFile = specLibFile;
		initialize();
	}
	
	private void initialize(){
		this.initializeIonIndexTable();
		this.table = initializeTable();
		this.initializeErrorModel();
	}
	
	private LookUpTable initializeTable(){
		LookUpTable table = new LookUpTable( 
			new int[] {lengthInterval.length, this.commonMaxCharge,  this.commonMaxCharge, this.sector.length, 
					this.ionsType.length+1, rankInterval.length+1, this.massToleranceInterval.length}); //extra slot at the beginning for noise model
		return table;
	}
	
	private void initializeErrorModel(){
		this.errorModel = new LookUpTable(
			new int[] {	this.rankInterval.length, this.massErrorInterval.length});    // assumed error model indenpendent of ranks
	}
	
	private void initializeIonIndexTable(){
		this.ionIndex = new HashMap<String, Integer>();
	//	this.ionsType = new String[Mass.standardIonsType.length+1];
		for(int i = 0; i < ionsType.length; i++){
			ionIndex.put(ionsType[i], new Integer(i));
			//System.out.println("storing ion: " + ionsType[i]);
		}
	}
	
	public int getIonIndex(LabelledPeak lp){
		if(!this.ionIndex.containsKey(lp.getType())){
			throw new IllegalArgumentException("Invalide ion type " + lp.getType());
		}
		return ionIndex.get(lp.getType()).intValue();
	}
	
	
	public double getValue(int[] index){
		return this.table.get(index);
	}
	//scoring table dimesnion
	//peaktype: common or linked
	//peakcharge
	//pepCharge
	//ion type
	private int[] getIndex(LabelledPeak lp, Peak realPeak){
		int rankIndex;
		int errorIndex;
		int sector;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
			sector = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
			errorIndex = ArrayUtils.getIntervalIndex(
					lp.getMass() - realPeak.getMass(), this.massToleranceInterval); 
			sector = realPeak.getSector();
		}
		if(lp == null){
			return new int[]{0,0,0,0,rankIndex};
		}
		Peptide p = lp.getPep();
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		///System.out.println("peptide is: " + p.getPeptide() + " index : " + peptideLength);
		int peptideCharge = lp.getPep().getCharge()-1;
		int peakCharge = lp.getCharge()-1;
		int ionIndex = getIonIndex(lp)+1;
		if(peptideCharge > 3){
			peptideCharge = 3;
		}
		if(peakCharge > 3){
			peakCharge = 3;
		}
		//System.out.println(peptideLength+"\t"+peptideCharge+"\t"+peakCharge+"\t"+ionIndex+"\t"+rankIndex+"\t"+errorIndex);
		return new int[]{peptideLength, peptideCharge, peakCharge, sector, ionIndex, rankIndex, errorIndex};
	}
	
	public int[] getErrorIndex(LabelledPeak lp, Peak realPeak){
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

	private int[] getNoiseIndex(Peptide p, Peak realPeak){
		int rankIndex;
		int errorIndex;
		int sector;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
			sector = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
			sector = realPeak.getSector();
		}
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		int peptideCharge = p.getCharge()-1;
		int peakCharge = 0;
		int ionIndex = 0;
		return new int[]{peptideLength, peptideCharge, peakCharge, sector, ionIndex, rankIndex, 0};
	}
	
	private int[] getNoiseIndex(LabelledPeak lp, Peak realPeak){
		int rankIndex;
		int errorIndex;
		int sector;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
			sector = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
			errorIndex = ArrayUtils.getIntervalIndex(
					lp.getMass() - realPeak.getMass(), this.massToleranceInterval); 
			sector = realPeak.getSector();
			
		}
		if(lp == null){
			return new int[]{0,0,0,sector, 0,rankIndex};
		}
		Peptide p = lp.getPep();
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		///System.out.println("peptide is: " + p.getPeptide() + " index : " + peptideLength);
		int peptideCharge = lp.getPep().getCharge()-1;
		int peakCharge = lp.getCharge()-1;
		int ionIndex = 0;
		return new int[]{peptideLength, peptideCharge, peakCharge, ionIndex, rankIndex, errorIndex};
	}
	
	public void getIonsCount(){
		String tripletFile ="..\\mixture_linked\\triplet_xquest.txt";
		LookUpTable totalCount = initializeTable();
		List<Spectrum> list = this.annotatedSet.getSpectrumList();
		for(int i = 0; i < list.size(); i++){
			Spectrum s = (Spectrum)list.get(i);
			s.windowFilterPeaks(12, 25);
			//s.computePeakSector();
			s.computePeakRank();
			//System.out.println("peptide is: " + s.peptide);
			//LabelledPeakFactory.setPoolMode(true);
			LabelledPeakFactory.resetFactory();
//			DecoySpectrumGenerator d = new DecoySpectrumGenerator();
			TheoreticalSpectrum t = new TheoreticalSpectrum(new Peptide(s.peptide+"."+s.charge), Mass.standardPrefixes, Mass.standardSuffixes);//+"."+s.charge));//, Mass.standardPrefixes, Mass.standardSuffixes);
			//System.out.println("pep " + s.peptide);
//			if(s.peptide.contains("+")|| s.peptide.contains("-")){
//				continue;
//			}
//			TheoreticalSpectrum t = new TheoreticalSpectrum(new Peptide(d.shuffle(s.peptide)+"."+s.charge), Mass.standardPrefixes, Mass.standardSuffixes);//+"."+s.charge));//, Mass.standardPrefixes, Mass.standardSuffixes);
//			TheoreticalSpectrum t = new TheoreticalSpectrum(new Peptide(s.peptide), Mass.prefixes_plus_noises, Mass.standardSuffixes);
			SimpleMatchingGraph matchingG = t.getMatchGraph(s, 0.5);
			this.getIonsCount(matchingG, totalCount);
		}
		this.normalizeCount();
		//this.smoothNoiseCount();
		this.normalizeErrorModel();
		this.printIonTable();
		this.printErrorTable();
		this.annotatedSet = null; //after training should release the training data and not to eat up memory
	}
	
	private Iterator<Spectrum> getAnnotatedIterator(){
		if(this.annotatedSet != null){
			return this.annotatedSet.getSpectrumList().iterator();
		}else{
			return new LargeSpectrumLibIterator(this.annnotateSetFile);
		}
	}
	
	public void getMixtureIonCount(){
		LookUpTable totalCount = initializeTable();
		int count = 0;
		for(Iterator it = getAnnotatedIterator(); it.hasNext();){
			Spectrum s = (Spectrum)it.next();
			s.computePeakRank();
			String[] peps = s.peptide.split(" & " );
			//System.out.println("peptide is: " + s.peptide);
			TheoreticalSpectrum t = new TheoreticalSpectrum(peps[0], peps[1]);
			SimpleMatchingGraph matchingG = t.getMatchGraph(s, 0.5);
			this.getIonsCount(matchingG, totalCount);
			count++;
			if(count % 1000 == 0){
				System.out.println("Finish Analyzing " + count);
			}
			//return;
		}
		this.normalizeCount();
		printIonTable();
	}
	
	public void getLinkedIonsCount(){
		String tripletFile ="..\\mixture_linked\\triplet_xquest.txt";
		LookUpTable totalCount = initializeTable();
		//String tripletFile =".\\mixture_linked\\triplet_selectedsubset.txt";
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); //skip header
			String[] tokens;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				System.out.println("line is: " + currentLine);
				tokens = currentLine.split("\t");
				System.out.println("id is : " + tokens[0]+".raw");
				m = this.annotatedSet.getSpectrumById(tokens[0]+".raw");
				System.out.println(m.peptide);
				m.computePeakRank();
				currentLine = bf.readLine();
				System.out.println("peptides are: " + tokens[1] + " & " + tokens[2]);
				//testMultipleLysPosition(tokens[1]+".2", tokens[2]+".2", m);
				Peptide p1 = new Peptide(tokens[1]+".2");
				Peptide p2 = new Peptide(tokens[2]+".2");
				p1.createDSPLinkerPTM(new int[]{Integer.parseInt(tokens[3])});
				p2.createDSPLinkerPTM(new int[]{ Integer.parseInt(tokens[4])});
				th = new TheoreticalSpectrum(p1, p2, (short)m.charge);
				SimpleMatchingGraph matchingG = th.getMatchGraph(m, 0.5);
				this.getIonsCount(matchingG, totalCount);
				//return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
		this.normalizeCount();
		printIonTable();
	}
	
	
	private void getIonsCount(SimpleMatchingGraph g, LookUpTable totalCount){
		Set vertices = g.vertexSet(2);
		Iterator it = vertices.iterator();
		Peak p, realPeak;
		LabelledPeak lp;
		Integer c, one = new Integer(1);
		Peptide pep = null;
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (LabelledPeak)p;
				pep = lp.getPep();
				Set<Peak> neighbors = g.getNeighborSet(lp);
				if(neighbors.size() == 0){
					int[] index = this.getIndex(lp, null);
					this.table.incrementIonCount(index);
				}
			}
		}
		
		it = g.vertexSet(1).iterator();
		while(it.hasNext()){
			p = (Peak)it.next();
			Set<Peak> neighbors = g.getNeighborSet(p);
			if(neighbors.size() == 0){				
				//int[] index = this.getIndex(null, p);
				int[] index = this.getNoiseIndex(pep, p);
				this.table.incrementIonCount(index);
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
				}
				for(Iterator<Peak> iter = neighbors.iterator(); iter.hasNext();){
					Peak neigh = iter.next();
					massDiff = Math.abs(neigh.getMass() - p.getMass());
					if(min == massDiff){
						int[] index = this.getIndex((LabelledPeak)neigh, p);
						this.table.incrementIonCount(index);
						int[] errorIndex = this.getErrorIndex((LabelledPeak)neigh, p);
						this.errorModel.incrementIonCount(errorIndex);
					}
				}
			}
		}
		//building noise model
		
	}
	
	//transform the ion counts to probabilities
	private void normalizeCount(){
		for(int length = 0; length < this.lengthInterval.length; length++){
			for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
				for(int peakCharge = 0; peakCharge < this.commonMaxCharge; peakCharge++){
					for(int sector = 0; sector < this.sector.length; sector++){
					for(int ionIndex = 0; ionIndex < this.ionsType.length+1; ionIndex++){
						double sum = 0.0;
						for(int rank = 0; rank < this.rankInterval.length+1; rank++){
							for(int noise = 0; noise < this.massToleranceInterval.length; noise++){
								int[] index = {length, pepCharge, peakCharge, sector, ionIndex, rank, noise};
								sum += this.table.get(index);
							}
						}
						for(int rank = 0; rank < this.rankInterval.length+1; rank++){
							for(int noise = 0; noise < this.massToleranceInterval.length; noise++){
								int[] index = {length, pepCharge, peakCharge, sector, ionIndex, rank, noise};
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
		//normalize for noises
		
		for(int length = 0; length < this.lengthInterval.length; length++){
			for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
				for(int sector = 0; sector < this.sector.length; sector++){
				double sum = 0.0;
				for(int rank = 0; rank < this.rankInterval.length+1; rank++){
					int[] index = {length, pepCharge, 0, sector, 0, rank, 0};
					sum += this.table.get(index);
				}
				sum /= (1-0.93);
				int[] noiseInd = {length, pepCharge, 0, sector, 0, 0,0};
				this.table.put(noiseInd, sum*0.95);
				for(int rank = 0; rank < this.rankInterval.length+1; rank++){
					int[] index = {length, pepCharge, 0,sector, 0, rank, 0};
					double count = this.table.get(index);
					this.table.put(index, count/sum);
				}
			  }
			}
		}
	}
	
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
	
	private void smoothNoiseCount(){
		for(int length = 0; length < this.lengthInterval.length-1; length++){
			for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
				for(int peakCharge = 0; peakCharge < this.commonMaxCharge; peakCharge++){
					for(int rank = 0; rank < this.rankInterval.length+1; rank++){							
						for(int noise = 0; noise < this.massToleranceInterval.length; noise++){
							double sum = 0.0;
							for(int ionIndex = this.ionsType.length - 2; ionIndex < this.ionsType.length+1; ionIndex++){
								int[] index = {length, pepCharge, peakCharge, ionIndex, rank, noise};
								if(rank == 0){
									//System.out.println("missing noise prob is: " + this.table.get(index));
								}
								sum += this.table.get(index);
							}
							int[] index = {length, pepCharge, peakCharge, 0, rank, noise};
							sum = sum / 3;
							if(rank == 0){
								///System.out.println("sum is: " + sum);
							}
							this.table.put(index, sum);	
						}
					}
				}
				int[] index = {length, pepCharge, 0, 0, 0, 0};
				//System.out.println("Average noise probability: "  + this.table.get(index));
			}
		}
	}
	
	private void printIonTable(){
		for(int length = 0; length < this.lengthInterval.length-1; length++){
			for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
				for(int peakCharge = 0; peakCharge < this.commonMaxCharge; peakCharge++){
					for(int sector = 0; sector < this.sector.length; sector++){
					for(int ionIndex = 0; ionIndex < this.ionsType.length; ionIndex++){
						String label;
						if(length == 0){
							label = "Short ";
						}else{
							label = "Long ";
						}
						double total = 0.0;
						StringBuffer header = new StringBuffer(label+ this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1) + " sector: " + sector);
						StringBuffer prob = new StringBuffer(header + " prob:\t");
						StringBuffer background = new StringBuffer(header + " noise:\t");
						for(int rank = 0; rank < this.rankInterval.length; rank++){
							for(int noise = 0; noise < this.massToleranceInterval.length-1; noise++){
								int[] index = {length, pepCharge, peakCharge, sector, ionIndex+1, rank+1, noise};
										prob.append(this.table.get(index) + "\t"); 
										background.append(this.table.get(new int[]{length, pepCharge, 0, sector, 0, rank+1, 0}) + "\t");
										total += this.table.get(index);
							}
						}
						System.out.println(prob.toString());
						System.out.println(background.toString());
						System.out.println(label+ this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1) + " sector: " + sector + 
								" total:  "      + total); 
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

	@Override
	public double compare(Peak p1, Peak p2) {
		LabelledPeak lp = (LabelledPeak)p1;
		if(p1 == null ){
			return 0;
		}else{
			//System.out.print("scoring: " + lp.getCharge() + "@" + lp.getCharge() + "@" + lp.getPep().getCharge());
			int[] index =  getIndex(lp, p2);
			int[] index2 = getNoiseIndex(lp.getPep(), p2); //assuem if we were to match the peak to noise
			int[] errorIndex = getErrorIndex(lp, p2);
//			int[] index2 = getNoiseIndex(lp, p2); //assuem if we were to match the peak to noise
			double score = this.table.get(index);
			double score2 = this.table.get(index2);
			if(p2 != null){
				double errorScore = this.errorModel.get(errorIndex);
				//System.out.print("peptide " +  lp +  "\t" + p2 + "\tscore: " + score + "\t" + score2 + "\t" + errorScore);
				//score *= errorScore;
				if(score == 0){
					score = 0.00001;
				}
				//score2 = score2 / (this.massErrorInterval.length-1);
				//System.out.println("\t" +  (1.0/(this.massErrorInterval.length-1)) + "\t" + Math.log(score/score2));
			}else{
				//System.out.println("peptide " +  lp +  "\t" + "-1.0\t-1.0\trank:\t-1.0" + "\tscore: " + score + "\t" + score2 + "\t" + 0.0 +"\t" +0.0 +"\t" + Math.log(score/score2));
			}
			//System.out.print("\tscore: " + score);
			//System.out.println("\tscore2: " + score2);
			if(Double.isNaN(score)){
				return 0;
			}
			if(score2 == 0){
				//System.out.println("warning noise score is zero");
				return 0;
				//score2 = 0.000001;
			}
			if(score < score2){
				//System.out.println("score is: " + score + " score2 is: " + score2);
				//return 0.0;
			}
			double ratio = Math.log(score/score2);
			int[] indexMissing = getIndex(lp, null);
			double scoreMissing = this.table.get(indexMissing);
			int[] indexMissing2 = getNoiseIndex(lp.getPep(), null);
			double scoreMissing2 = this.table.get(indexMissing2);
			double ratio2 = Math.log(scoreMissing/scoreMissing2);
			//System.out.println("ratio1:\t" + ratio + "\tratio2\t" + ratio2);
			return ratio;
//			if(ratio2 > ratio){
//				return ratio2; 
//			}else{
//				return ratio;
//			}
		}
	}
	
	public static List<Spectrum> generateSpectra(List<String> pepList, Spectrum linkedquery){
		List<Spectrum> specList = new ArrayList<Spectrum>();
		for(Iterator<String> it = pepList.iterator(); it.hasNext();){
			String pep = it.next();
			Peptide p = new Peptide(pep+".2");
			p.setPtmmasses(new double[]{LookUpSpectrumLibX.getLinkedOffSet(pep, linkedquery)});
			int pos = pep.indexOf('K');
			while(pos > 0){
				Peptide copy = new Peptide(p);
				p.setPos(new int[]{pos});
				TheoreticalSpectrum th = new TheoreticalSpectrum(p, linkedquery.charge);
				specList.add(th);
				pos = pep.indexOf('K', pos+1);
			}
		}
		return specList;
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
	
	public static RankBaseScoreLearner loadComparator(String file){
		try{
			BufferedInputStream bi = new BufferedInputStream(new FileInputStream(file));
			ObjectInputStream oi = new ObjectInputStream(bi);
		    Object o = oi.readObject();
		    return (RankBaseScoreLearner)o;
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	public static void testGetIonStat(){
		String file = "..\\mixture_linked\\test.mgf";
		SpectrumLib lib = new SpectrumLib(file, "MGF");
		//lib.removeModSpectra();
		lib.windowFilterPeaks(12, 25);
		lib.computeRank();
		RankBaseScoreLearner learner = new RankBaseScoreLearner(lib);
		learner.getIonsCount();
		learner.writeLibToFile("..\\mixture_linked\\Lib_disulfide_dangle.o");
	}
	
	
	public static void testScorer(){
		String file = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib = new SpectrumLib(file, "MSP");
		String pepfile = "..\\mixture_linked\\Ecoli_allpeptides_plusLinkedpeptides.txt";
		lib.removeModSpectra();
		lib.computeRank();
		//lib.windowFilterPeaks(6, 25);
		RankBaseScoreLearner learner = new RankBaseScoreLearner(lib);
		learner.getIonsCount();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(pepfile);
		List<Spectrum> testList = SpectrumUtil.getRandomSpectrum(lib, 500, 3);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		SpectrumIonRankLearner learner2 = new SpectrumIonRankLearner(lib);
		PeakComparator peakscorer2 = learner2.createComparatorSet();
		SimpleProbabilisticScorer scorer2 = new SimpleProbabilisticScorer(peakscorer2);
		for(int i = 0; i < testList.size(); i++){
			Spectrum query = testList.get(i);
			query.windowFilterPeaks(6, 25);
			SpectrumLib cand = factory.createCandidateSpectrumLibX(query, 3, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cand.getSpectrumList(), scorer, scorer);
			SpectrumLibSearcher searcher2 = new SpectrumLibSearcher(cand.getSpectrumList(), scorer2, scorer2);
			searcher.topSpectra(query, 5);
			searcher2.topSpectra(query, 5);
		}
	}
	
	public static void testScoreFilter(){
		String libfile = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib = new SpectrumLib(libfile, "MSP");
		lib.removeModSpectra();
		LinkedPeakScoreLearner learner = new LinkedPeakScoreLearner(lib);
		learner.getIonsCount();
		String file = "..\\mixture_linked\\Ecoli_allpeptides_plusLinkedpeptides.txt";
		//String file = "..\\mixture_linked\\tempLinkedPeptides.txt";
		List<Spectrum> specList = lib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		LookUpSpectrumLibX lookup = new LookUpSpectrumLibX(factory.peptides);
		factory = null;
		System.out.println("Done indexing peptides");
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		scorer.includeNoise = false;
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			List<Peak> queryPeaks = s.getTopPeaks(20);
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<String> candidates = lookup.getSpectrumByPeaks(queryPeaks, 2, s);
			String[] peptides = s.peptide.split(" & ");
			int passedFilter = lookup.checkPassFilter(peptides[0], peptides[1].substring(0, peptides[1].length()-2), candidates);
			System.out.println("Query: "  + s.spectrumName + " After filter one we have: " + candidates.size() + " candidates ");
			System.out.println("After filter correct peptide is retained?: " + passedFilter);
			List<Spectrum> candidateSpectrum = generateSpectra(candidates, s);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidateSpectrum, scorer);
			int[] ranks = searcher.linkedRanks(s);
			System.out.println("target peptides ranks " + ranks[0] + "\t" + ranks[1]);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		
	}
	
	public static void main(String[] args){
		testGetIonStat();
		//testScoreFilter();
		//testScorer();
	}

	public int[] getRankInterval() {
		return rankInterval;
	}

	public void setRankInterval(int[] rankInterval) {
		this.rankInterval = rankInterval;
	}

	public int[] getLengthInterval() {
		return lengthInterval;
	}

	public void setLengthInterval(int[] lengthInterval) {
		this.lengthInterval = lengthInterval;
	}
}
