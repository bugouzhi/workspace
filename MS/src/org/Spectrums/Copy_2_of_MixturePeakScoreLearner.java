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
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Copy_2_of_MixturePeakScoreLearner implements PeakComparator, Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1956042383355436987L;
	private String[] ionsType = Mass.standardIonsType;
	private static int MAXRANK = 2000;
	private static int MAXLENGTH = 150;
	private static int PEPCOUNT = 2;
	private int[] rankInterval = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30, 35, 
			40, 45, 50, 55,60,65, 70,75,80, 85, 90,100,110, 120, 130, 140, 150,160,170,180,190,200,210,220,230,240,250, MAXRANK};
	private int[] lengthInterval = {1,12,MAXLENGTH};
	//to define an convention, here mass error means theoretical mass - actual mass of a peak
	private double[] massToleranceInterval = {-0.5,0.5}; //need to setup this automatically in the future
	private double[] massErrorInterval = {-0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55};
	private HashMap<String, Integer> ionIndex;
	private int commonMinCharge = 1;
	private int commonMaxCharge = 5;
	private int linkedMinCharge = 2;
	private int linkedMaxCharge = 8;
	private Map<String, Double> priorProbs;
	private LookUpTable table;
	private LookUpTable errorModel;
	private SpectrumLib annotatedSet;
	private String annnotateSetFile;
	
	public Copy_2_of_MixturePeakScoreLearner(String file){
		this.annnotateSetFile = file;
		initialize();
	}
	
	
	private void initialize(){
		this.initializeIonIndexTable();
		this.table = initializeTable();
		this.initializeErrorModel();
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
		
	private LookUpTable initializeTable(){
		LookUpTable table = new LookUpTable( 
			new int[] {this.PEPCOUNT, this.linkedMaxCharge, lengthInterval.length, this.commonMaxCharge,  this.commonMaxCharge, 
					this.ionsType.length+1, rankInterval.length+1, this.massToleranceInterval.length}); //extra slot at the beginning for noise model
		return table;
	}
	
	@Override
	public double compare(Peak p1, Peak p2) {
		MixturePeak lp = (MixturePeak)p1;
		if(p1 == null ){
			return 0;
		}else{
			//System.out.println("scoring: " + lp.getCharge() + "@" + lp.getCharge() + "@" + lp.getPep().getCharge());
			int[] index =  getIndex(lp, p2);
			int[] index2 = getNoiseIndex(lp.getPep(), p2, lp.getParent().charge, lp.getPeptideIndex()); //assuem if we were to match the peak to noise
			int[] errorIndex = getErrorIndex(lp, p2);
//			int[] index2 = getNoiseIndex(lp, p2); //assuem if we were to match the peak to noise
			double score = this.table.get(index);
			double score2 = this.table.get(index2);
			if(p2 != null){
				//System.out.println("peptide " + lp.getPeptideIndex()  +  "\t" +  lp +  "\trank: " + p2.getRank() + " score: " + score + "\t" + score2 + "\t" + Math.log(score/score2));
				double errorScore = this.errorModel.get(errorIndex);
				score *= errorScore;
				if(score == 0){
					score = 0.00001;
				}
				score2 = score2 / (this.massErrorInterval.length-1);
			}
//			System.out.println("score: " + score);
//			System.out.println("score2: " + score2);
			if(Double.isNaN(score)){
				return 0;
			}
			if(score2 == 0){
				System.out.println("warning noise score is zero");
				System.out.println("index is: " + lp.getParent().charge + "\t" + lp.getPeptideIndex() + "\t"
						+ lp.getPep() + "\t" + lp.getType() + "\t" + lp.getCharge());
				score2 = 0.00001;
			}
			if(score < score2){
				//System.out.println("score is: " + score + " score2 is: " + score2);
				//return 0.0;
			}
			return Math.log(score/score2);
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
	//peptide index
	//combined charge
	//length
	//peakcharge
	//pepCharge
	//ion type
	//mass error
	
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
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		///System.out.println("peptide is: " + p.getPeptide() + " index : " + peptideLength);
		int peptideCharge = lp.getPep().getCharge()-1;
		int peakCharge = lp.getCharge()-1;
		int ionIndex = getIonIndex(lp)+1;
		return new int[]{lp.getPeptideIndex(), lp.getParent().charge, peptideLength, peptideCharge, peakCharge, ionIndex, rankIndex, errorIndex};
	}
	
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

	private int[] getNoiseIndex(Peptide p, Peak realPeak, int combineCharge, int peptideIndex){
		int rankIndex;
		int errorIndex;
		if(realPeak == null){
			rankIndex = 0;
			errorIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
		}
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.lengthInterval);
		int peptideCharge = p.getCharge()-1;
		int peakCharge = 0;
		int ionIndex = 0;
		return new int[]{peptideIndex, combineCharge, peptideLength, peptideCharge, peakCharge, ionIndex, rankIndex, 0};
	}
	
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
	
	public void getMixtureIonCount(){
		LookUpTable totalCount = initializeTable();
		int count = 0;
		for(Iterator it = getAnnotatedIterator(); it.hasNext();){
			Spectrum s = (Spectrum)it.next();
			if(s.modMass > 0){
				continue;
			}
			s.windowFilterPeaks(12, 50);
			s.computePeakRank();
			String[] peps = s.peptide.split(" & " );
			//System.out.println("peptide is: " + s.peptide);
			TheoreticalSpectrum t = new TheoreticalSpectrum(peps[0], peps[1]);
			SimpleMatchingGraph matchingG = t.getMatchGraph(s, 0.5);
			this.getIonsCount(matchingG, peps[0], peps[1]);
			count++;
			if(count % 1000 == 0){
				System.out.println("Finish Analyzing " + count);
			}
			//return;
//			if(count == 10000){
//				break;
//			}
		}
		this.normalizeCount();
		this.normalizeErrorModel();
		printIonTable();
		printErrorTable();
	}
	
	private void normalizeCount(){
		for(int peptide = 0; peptide < this.PEPCOUNT; peptide++){
			for(int combineCharge = 0; combineCharge < this.linkedMaxCharge; combineCharge++){
				for(int length = 0; length < this.lengthInterval.length; length++){
					for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
						for(int peakCharge = 0; peakCharge < this.commonMaxCharge; peakCharge++){
							for(int ionIndex = 0; ionIndex < this.ionsType.length+1; ionIndex++){
								double sum = 0.0;
								for(int rank = 0; rank < this.rankInterval.length+1; rank++){
									for(int noise = 0; noise < this.massToleranceInterval.length; noise++){
										int[] index = {peptide, combineCharge, length, pepCharge, peakCharge, ionIndex, rank, noise};
										sum += this.table.get(index);
									}
								}
								for(int rank = 0; rank < this.rankInterval.length+1; rank++){
									int binWidth = 1;
									if(rank < this.rankInterval.length && rank > 0){
										binWidth = this.rankInterval[rank]-this.rankInterval[rank-1];
									}	
									for(int noise = 0; noise < this.massToleranceInterval.length; noise++){
										int[] index = {peptide, combineCharge, length, pepCharge, peakCharge, ionIndex, rank, noise};
										double count = this.table.get(index);
										if(count > sum) System.out.println("waring sum is smaller: " + count + " out of " + sum);
										this.table.put(index, count/binWidth/sum);
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
			for(int combineCharge = 0; combineCharge < this.linkedMaxCharge; combineCharge++){
				for(int length = 0; length < this.lengthInterval.length; length++){
					for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
						double sum = 0.0;
						for(int rank = 0; rank < this.rankInterval.length+1; rank++){
							int[] index = {peptide, combineCharge, length, pepCharge, 0, 0, rank, 0};
							sum += this.table.get(index);
						}
						sum /= (1-0.93);
						int[] noiseInd = {peptide, combineCharge, length, pepCharge, 0, 0, 0,0};
						this.table.put(noiseInd, sum*0.95);
						for(int rank = 0; rank < this.rankInterval.length+1; rank++){
							int binWidth=1;
							if(rank < this.rankInterval.length && rank > 0){
								binWidth = this.rankInterval[rank]-this.rankInterval[rank-1];
							}
							int[] index = {peptide, combineCharge, length, pepCharge, 0, 0, rank, 0};
							double count = this.table.get(index);
							this.table.put(index, count/binWidth/sum);
						}
					}
				}
			}
		}
	}
	
	private void getIonsCount(SimpleMatchingGraph g, String pep1, String pep2){
		Set vertices = g.vertexSet(2);
		Iterator it = vertices.iterator();
		Peak p, realPeak;
		MixturePeak lp;
		Integer c, one = new Integer(1);
		Peptide pep = null;
		String[] peps1 = pep1.split("\\.");
		String[] peps2 = pep2.split("\\.");
		Peptide peptide1 = new Peptide(pep1);
		Peptide peptide2 = new Peptide(pep2);
		int combineCharge = Integer.parseInt(peps1[1]) 
			+ Integer.parseInt(peps2[1]);
		int peptideIndex = 0;
		
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (MixturePeak)p;
				pep = lp.getPep();
				Set<Peak> neighbors = g.getNeighborSet(lp);
				if(lp.getPep().getPeptide().equals(peps1[0])){
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
			Set<Peak> neighbors = g.getNeighborSet(p);
			if(neighbors.size() == 0){				
				//int[] index = this.getIndex(null, p);
				int[] index = this.getNoiseIndex(peptide1, p, combineCharge, 0);
				int[] index2 = this.getNoiseIndex(peptide2, p, combineCharge, 1);
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
				}
				if(((LabelledPeak)closest).getPep().getPeptide().equals(peps1[0])){
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
		//building noise model
		
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
	
	private void printIonTable(){
		for(int peptide = 0; peptide < this.PEPCOUNT; peptide++){
			for(int combineCharge = 0; combineCharge < this.linkedMaxCharge; combineCharge++){
				for(int length = 0; length < this.lengthInterval.length-1; length++){
					for(int pepCharge = 0; pepCharge < this.commonMaxCharge; pepCharge++){
						for(int peakCharge = 0; peakCharge < this.commonMaxCharge; peakCharge++){
							for(int ionIndex = 0; ionIndex < this.ionsType.length; ionIndex++){
								for(int rank = 0; rank < this.rankInterval.length; rank++){
									for(int noise = 0; noise < this.massToleranceInterval.length-1; noise++){
										int[] index = {peptide, combineCharge, length, pepCharge, peakCharge, ionIndex+1, rank+1, noise};
										if(length == 0){
											System.out.print("Short(<14) : ");
										}else{
											System.out.print("Long(>14) : ");
										}
										String label;
										if(peptide == 0){
											label = "Peptide1 ("+ combineCharge + ") ";
										}else{
											label = "Peptide2 (" + combineCharge + ") ";
										}
										System.out.println(label+ this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1) + 
												" rank " + this.rankInterval[rank] +  " error: " + noise + ": "    
												+ this.table.get(index) 
												+ " noise: " + this.table.get(new int[]{peptide, combineCharge, length, pepCharge, 0, 0, rank+1, 0}));
									}
								}
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
	
	public static Copy_2_of_MixturePeakScoreLearner loadComparator(String file){
		try{
			BufferedInputStream bi = new BufferedInputStream(new FileInputStream(file));
			ObjectInputStream oi = new ObjectInputStream(bi);
		    Object o = oi.readObject();
		    return (Copy_2_of_MixturePeakScoreLearner)o;
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	public static void testLoadComparator(){
		String trainFile = "..\\mixture_linked\\mixtures100000_alpha1.0.mgf";
		String outfile = "..\\mixture_linked\\mixtures_alpha_models.o";
		Copy_2_of_MixturePeakScoreLearner peakscorer = new Copy_2_of_MixturePeakScoreLearner(trainFile); //scorer
		peakscorer.getMixtureIonCount();
		peakscorer.writeLibToFile(outfile);
		Copy_2_of_MixturePeakScoreLearner peakscorer2 = loadComparator(outfile);
		System.out.println(peakscorer2);
		
	}
	public static void testMixtureScoring(){
		String outfile = "..\\mixture_linked\\mixtures_generic_models.o";
		String outfile2 = "..\\mixture_linked\\mixtures_alpha0.3_models.o";
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib.txt";
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		String mixturefile = "..\\mixture_linked\\yeast_mixture.name"; 
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		
//		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
//		peakscorer2.getIonsCount();
//		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);

		Copy_2_of_MixturePeakScoreLearner peakscorer2 = loadComparator(outfile); //scorer
		MixtureSpectrumScorer scorer2 = new MixtureSpectrumScorer(peakscorer2);

		Copy_2_of_MixturePeakScoreLearner peakscorer3 = loadComparator(outfile2); //scorer
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(peakscorer3);
		
		//SpectrumLib mixlib = lib1.createMix(mixturefile, 1000, 1, 0.3, 0.0001, 1, 2, false);
		SpectrumLib mixlib = lib1.createRandomMix(1000, 1, 0.3, 0.0001, 1, 2, false);
		lib1 = null;
		List<Spectrum>  specList = mixlib.getSpectrumList();
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			String[] peptides = s.peptide.split(" & ");
			TheoreticalSpectrum t = new TheoreticalSpectrum(peptides[0], peptides[1]);
			double score1 = scorer.compare(t, s);
			double score2 = scorer2.compare(t, s);
			System.out.println(s.peptide + " has PSM score: " + score1 + " and " + score2);
			
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMixtureModel(){
		String trainFile = "..\\mixture_linked\\mixtures_generic_models.o";
		String training = "..\\mixture_linked\\mixtures100000_alpha0.5.mgf";
		//String training = "..\\mixture_linked\\yeast_mixture_alpha_high.mgf";
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib_plusDecoy.txt";
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		String mixturefile = "..\\mixture_linked\\yeast_mixture.name"; 
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);
		
		//MixturePeakScoreLearner peakscorer3 = loadComparator(trainFile); //scorer
		Copy_2_of_MixturePeakScoreLearner peakscorer3 = new Copy_2_of_MixturePeakScoreLearner(training);
		peakscorer3.getMixtureIonCount();
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(peakscorer3);
		
		//SpectrumLib mixlib = lib1.createMix(mixturefile, 100, 1, 0.5, 0.0001, 1, 2, false);
		SpectrumLib mixlib = lib2.createRandomMix(100, 0.5,  0.0001, 1.0, 2, false);
		List<Spectrum>  specList = mixlib.getSpectrumList();
//		List<Spectrum>  specList = new ArrayList<Spectrum>();
//		for(int i = 0; i < 5000; i++){
//			specList.add(lib1.getRandomSpectrum());
//		}
		lib1 = null;
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		factory.setMinCharge(2);
		factory.setMaxCharge(3);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(10, 50);
			s.computePeakRank();
			//s.peptide = s.peptide + " & " + s.peptide;
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), filter, scorer);
			searcher.setSingleScorer(scorer1);
			//int[] ranks = searcher.ranks(s);
			//System.out.println("Spectrum: " + s.peptide + " ranks " +  ranks[0] + " & " + ranks[1]+ " in spectrumLib of size: " + lib.getSpectrumList().size());
			searcher.bestCandidates(s, 10);
			//searcher.bestSpectrum(s);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	
	public static void testMixtureModelExperimental(){
		String trainFile = "..\\mixture_linked\\mixtures_generic_models.o";
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_122007p_yeast_digest1.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib_plusDecoy.txt";
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		System.out.println("Starting training");
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);
		//lib1 = null;
		Copy_2_of_MixturePeakScoreLearner peakscorer3 = loadComparator(trainFile); //scorer
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(peakscorer3);
		System.out.println("Done training");

//		List<Spectrum>  specList = new ArrayList<Spectrum>();
		LargeSpectrumLibIterator<Spectrum> iterator = new LargeSpectrumLibIterator<Spectrum>(spectrumFile2);
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(; iterator.hasNext();){
			Spectrum s = iterator.next();
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			//s.peptide = s.peptide + " & " + s.peptide;
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), filter, scorer);
			searcher.setSingleScorer(scorer1);
			//int[] ranks = searcher.ranks(s);
			//System.out.println("Spectrum: " + s.peptide + " ranks " +  ranks[0] + " & " + ranks[1]+ " in spectrumLib of size: " + lib.getSpectrumList().size());
			searcher.bestCandidates(s, 10);
			//searcher.bestSpectrum(s);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMixtureModelSubset(){
		String trainFile = "..\\mixture_linked\\mixtures_alpha_generic.o";
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
		String querySpectrum = "..\\mixture_linked\\mixtureQuery.txt";
		String training = "..\\mixture_linked\\mixtures100000_alpha0.5.mgf";
		//String training = "..\\mixture_linked\\yeast_mixture_alpha_high.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib_plusDecoy.txt";
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		System.out.println("Starting training");
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);
		//lib1 = null;
//		MixturePeakScoreLearner peakscorer3 = loadComparator(trainFile); //scorer
//		MixturePeakScoreLearner peakscorer3 = new MixturePeakScoreLearner(training);
//		peakscorer3.getMixtureIonCount();
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(peakscorer2);
		System.out.println("Done training");
		SpectrumLib mixLib = new SpectrumLib(spectrumFile2, "MGF");
		Iterator<Spectrum> iterator = SpectrumUtil.getSpectra(querySpectrum, mixLib).iterator();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		factory.setMinCharge(2);
		factory.setMaxCharge(3);
		
		long start = (new GregorianCalendar()).getTimeInMillis();
		mixLib = null;
		for(; iterator.hasNext();){
			Spectrum s = iterator.next();
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(12, 25);
			//s.windowFilterAndRank(7, 50, 150);
			s.computePeakRank();
			//s.peptide = s.peptide + " & " + s.peptide;
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), filter, scorer);
			searcher.setSingleScorer(scorer1);
			//int[] ranks = searcher.ranks(s);
			//System.out.println("Spectrum: " + s.peptide + " ranks " +  ranks[0] + " & " + ranks[1]+ " in spectrumLib of size: " + lib.getSpectrumList().size());
			searcher.bestCandidates(s, 10);
			//searcher.bestSpectrum(s);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	public static void main(String[] args){
		//testLoadComparator();
		//testMixtureScoring();
		testMixtureModel();
		//testMixtureModelExperimental();
		//testMixtureModelSubset();
	}

}
