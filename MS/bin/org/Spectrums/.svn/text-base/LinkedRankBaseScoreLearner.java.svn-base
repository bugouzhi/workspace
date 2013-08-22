package org.Spectrums;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.jgrapht.graph.SimpleGraph;

public class LinkedRankBaseScoreLearner implements PeakComparator{
	private String[] ionsType = Mass.standardIonsType;
	private static int MAXRANK = 2000;
	private int[] rankInterval = {1, 3, 10, 15, 20, 30, 60, 100, MAXRANK};
	private HashMap<String, Integer> ionIndex;
	private int commonMinCharge = 1;
	private int commonMaxCharge = 4;
	private int linkedMinCharge = 2;
	private int linkedMaxCharge = 6;
	private LookUpTable table;
	private SpectrumLib annotatedSet;
	public LinkedRankBaseScoreLearner(SpectrumLib lib){
		this.annotatedSet = lib;
		this.initializeIonIndexTable();
		this.table = initializeTable();
	}
	
	private LookUpTable initializeTable(){
		LookUpTable table = new LookUpTable( 
			new int[] {2, linkedMaxCharge,  linkedMaxCharge, this.ionsType.length+1, rankInterval.length+1}); //extra slot at the beginning for noise model
		return table;
	}
	
	private void initializeIonIndexTable(){
		this.ionIndex = new HashMap<String, Integer>();
	//	this.ionsType = new String[Mass.standardIonsType.length+1];
		for(int i = 0; i < ionsType.length; i++){
			ionIndex.put(ionsType[i], new Integer(i));
			//System.out.println("storing ion: " + ionsType[i]);
		}
	}
	
	private int getIonIndex(LabelledPeak lp){
		if(!this.ionIndex.containsKey(lp.getType())){
			throw new IllegalArgumentException("Invalide ion type " + lp.getType());
		}
		return ionIndex.get(lp.getType()).intValue();
	}
	
	//scoring table dimesnion
	//peaktype: common or linked
	//peakcharge
	//pepCharge
	//ion type
	private int[] getIndex(LabelledPeak lp, Peak realPeak){
		int peakType;
		int rankIndex;
		if(realPeak == null){
			rankIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.rankInterval)+1;
		}
		if(lp == null){
			return new int[]{0,0,0,0,rankIndex};
		}
		Peptide p = lp.getPep();
		if(TheoreticalSpectrum.isLinkedPeak(p, lp)){
			peakType = 1;
		}else{
			peakType = 0;
		}
		int peptideCharge = lp.getPep().getCharge()-1;
		int peakCharge = lp.getCharge()-1;
		int ionIndex = getIonIndex(lp)+1;
		return new int[]{peakType, peptideCharge, peakCharge, ionIndex, rankIndex};
	}
	
	public void getIonsCount(){
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
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (LabelledPeak)p;
				Set<Peak> neighbors = g.getNeighborSet(lp);
				if(neighbors.size() > 0){
					realPeak = neighbors.iterator().next();
				}else{
					realPeak = null;
				}
				int[] index = this.getIndex(lp, realPeak);
				this.table.incrementIonCount(index);
			}
		}
		
		it = g.vertexSet(1).iterator();
		while(it.hasNext()){
			p = (Peak)it.next();
			Set<Peak> neighbors = g.getNeighborSet(p);
			if(neighbors.size() == 0){
				int[] index = this.getIndex(null, p);
				this.table.incrementIonCount(index);
			}
		}
		//building noise model
		
	}
	
	//transform the ion counts to probabilities
	private void normalizeCount(){
		for(int type = 0; type < 2; type++){
			for(int pepCharge = 2; pepCharge < 6; pepCharge++){
				for(int peakCharge = 0; peakCharge < 6; peakCharge++){
					for(int ionIndex = 0; ionIndex < this.ionsType.length; ionIndex++){
						double sum = 0.0;
						for(int rank = 0; rank < this.rankInterval.length; rank++){
							int[] index = {type, pepCharge, peakCharge, ionIndex, rank};
							sum += this.table.get(index);
						}
						for(int rank = 0; rank < this.rankInterval.length; rank++){
							int[] index = {type, pepCharge, peakCharge, ionIndex, rank};
							double count = this.table.get(index);
							this.table.put(index, count/sum);
						}

					}
				}
			}
		}
		//normalize for noises
		double sum = 0.0;
		for(int rank = 0; rank < this.rankInterval.length; rank++){
			int[] index = {0, 0, 0, 0, rank};
			sum += this.table.get(index);
		}
		sum /= (1-0.95);
		for(int rank = 0; rank < this.rankInterval.length; rank++){
			int[] index = {0, 0, 0, 0, rank};
			double count = this.table.get(index);
			this.table.put(index, count/sum);
		}
	}
	
	private void printIonTable(){
		for(int type = 0; type < 2; type++){
			for(int pepCharge = 2; pepCharge < 6; pepCharge++){
				for(int peakCharge = 0; peakCharge < 6; peakCharge++){
					for(int ionIndex = 0; ionIndex < this.ionsType.length; ionIndex++){
						for(int rank = 0; rank < this.rankInterval.length; rank++){
							int[] index = {type, pepCharge, peakCharge, ionIndex+1, rank+1};
							if(type == 0){
								System.out.print("unlinked: ");
							}else{
								System.out.print("linked: ");
							}
							System.out.println(this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1) + " rank " + (rank+1) +  ": "  
									+ this.table.get(index) + " noise: " + this.table.get(new int[]{0,0,0,0, rank+1}));
						}
					}
				}
			}
		}
	}
	

	@Override
	public double compare(Peak p1, Peak p2) {
		LabelledPeak lp = (LabelledPeak)p1;
		if(p1 == null || p2 == null){
			return 0;
		}else{
			//System.out.println("scoring: " + lp.getCharge() + "@" + lp.getCharge() + "@" + lp.getPep().getCharge());
			int[] index =  getIndex(lp, p2);
			double score = this.table.get(index);
			if(Double.isNaN(score)){
				return 0;
			}
			if(score < 0.05){
				return 0;
			}
			//System.out.println("matching score is: " + score);
			//return 1.0;
			return Math.log(this.table.get(index)/0.05);
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
	
	public static void testGetIonStat(){
		String file = "..\\mixture_linked\\spectrums_raw.mgf";
		SpectrumLib lib = new SpectrumLib(file, "MGF");
		LinkedRankBaseScoreLearner learner = new LinkedRankBaseScoreLearner(lib);
		learner.getIonsCount();
		//learner.getNoisesStatistics();
	}
	
	public static void testScoreFilter(){
		String libfile = "..\\mixture_linked\\linked_peptide_spectra2.mgf";
		SpectrumLib lib = new SpectrumLib(libfile, "MGF");
		String libfile2 = "..\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib2 = new SpectrumLib(libfile2, "MSP");
		lib2.removeModSpectra();
		lib2.computeRank();
		LPeakRankBaseScorer learner = new LPeakRankBaseScorer(lib2);
		String file = "..\\mixture_linked\\Ecoli_allpeptides_plusLinkedpeptides.txt";
		//String file = "..\\mixture_linked\\tempLinkedPeptides.txt";
		List<Spectrum> specList = lib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		LookUpSpectrumLibX lookup = new LookUpSpectrumLibX(factory.peptides, 0.5);
		factory = null;
		System.out.println("Done indexing peptides");
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		scorer.includeNoise = false;
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			//s.windowFilterPeaks(10, 25);
			s.computePeakRank();
			List<Peak> queryPeaks = s.getTopPeaks(20);
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<String> candidates = lookup.getSpectrumByPeaks(queryPeaks, 3, s);
			String[] peptides = s.peptide.split(" & ");
			//int passedFilter = lookup.checkPassFilter(peptides[0], peptides[1].substring(0, peptides[1].length()-2), candidates);
			int passedFilter = lookup.checkPassFilter(peptides[0], peptides[1], candidates);
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
		//testGetIonStat();
		testScoreFilter();
	}

}
