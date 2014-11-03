package org.Spectrums;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.jgrapht.graph.SimpleGraph;

public class LinkedPeakScoreLearner implements PeakComparator{
	private String[] ionsType = Mass.standardIonsType;
	private HashMap<String, Integer> ionIndex;
	private int commonMinCharge = 1;
	private int commonMaxCharge = 4;
	private int linkedMinCharge = 2;
	private int linkedMaxCharge = 6;
	private LookUpTable table;
	private SpectrumLib annotatedSet;
	public LinkedPeakScoreLearner(SpectrumLib lib){
		this.annotatedSet = lib;
		this.initializeIonIndexTable();
		this.table = initializeTable();
	}
	
	private LookUpTable initializeTable(){
		LookUpTable table = new LookUpTable( 
			new int[] {2, linkedMaxCharge,  linkedMaxCharge, this.ionsType.length});
		return table;
	}
	
	private void initializeIonIndexTable(){
		this.ionIndex = new HashMap<String, Integer>();
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
		Peptide p = lp.getPep();
		if(TheoreticalSpectrum.isLinkedPeak(p, lp)){
			peakType = 1;
		}else{
			peakType = 0;
		}
		int peptideCharge = lp.getPep().getCharge()-1;
		int peakCharge = lp.getCharge()-1;
		int ionIndex = getIonIndex(lp);
		return new int[]{peakType, peptideCharge, peakCharge, ionIndex};
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
				SimpleGraph matchingG = th.matchSpectrum2(m, 0.5);
				this.getIonsCount(matchingG, totalCount);
				//return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}

		LookUpTable.mapOperator(this.table, totalCount, TableElementOperator.Divider.d);
		printIonTable();
	}
	
	
	private void getIonsCount(SimpleGraph g, LookUpTable totalCount){
		Set vertices = g.vertexSet();
		Iterator it = vertices.iterator();
		Peak p;
		LabelledPeak lp;
		Integer c, one = new Integer(1);
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (LabelledPeak)p;
				int[] index = getIndex(lp, lp);
				if(g.degreeOf(p) > 0){
					this.table.incrementIonCount(index);
				}
				totalCount.incrementIonCount(index);
			}
		}
	}
	
	private void printIonTable(){
		for(int type = 0; type < 2; type++){
			for(int pepCharge = 2; pepCharge < 6; pepCharge++){
				for(int peakCharge = 0; peakCharge < 6; peakCharge++){
					for(int ionIndex = 0; ionIndex < this.ionsType.length; ionIndex++){
						int[] index = {type, pepCharge, peakCharge, ionIndex};
						if(type == 0){
							System.out.print("unlinked: ");
						}else{
							System.out.print("linked: ");
						}
						System.out.println(this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1) + ": " + this.table.get(index));
				
					}
				}
			}
		}
	}
	

	@Override
	public double compare(Peak p1, Peak p2) {
		LabelledPeak lp = (LabelledPeak)p1;
		if(p1 == null){
			return 0;
		}if(p2 == null){
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
			return Math.log((1-this.table.get(index))/(1-0.05));
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
			while(pos >= 0){
				Peptide copy = new Peptide(p);
				copy.setPos(new int[]{pos+1});
//				TheoreticalSpectrum th = new TheoreticalSpectrum(copy, linkedquery.charge);
				TheoreticalSpectrum th = new LazyEvaluateLinkedSpectrum(copy, linkedquery.charge);
				//TheoreticalSpectrum th = new TheoreticalSpectrum(p);
				specList.add(th);
				pos = pep.indexOf('K', pos+1);
			}
		}
		return specList;
	}
	
	
	
	public static List<Peptide> generatePeptides(List<PeptideLite> pepList){
		List<Peptide> peptides = new ArrayList<Peptide>();
		for(Iterator<PeptideLite> it = pepList.iterator(); it.hasNext();){
			PeptideLite current = it.next();
			String pep = current.getPep();
			//System.out.println("peptide is: " + pep);
			Peptide p = new Peptide(pep+".1");
			peptides.add(p);
		}
		return peptides;
	}
	
	public static List<Spectrum> generateSpectra3(List<Peptide> pepList, Spectrum linkedquery){
		List<Spectrum> specList = new ArrayList<Spectrum>();
		for(Iterator<Peptide> it = pepList.iterator(); it.hasNext();){
			Peptide p = it.next();
			//System.out.println(p);
			TheoreticalSpectrum th = new LazyEvaluateLinkedSpectrum(p, linkedquery.charge);
			specList.add(th);
		}
		return specList;
	}
	
	public static List<Spectrum> generateSpectra2(List<PeptideLite> pepList, Spectrum linkedquery){
		List<Spectrum> specList = new ArrayList<Spectrum>();
		for(Iterator<PeptideLite> it = pepList.iterator(); it.hasNext();){
			PeptideLite current = it.next();
			String pep = current.getPep();
			Peptide p = new Peptide(pep+".2");
			p.setPtmmasses(new double[]{LookUpSpectrumLibX.getLinkedOffSet(pep, linkedquery)});
			int pos = pep.indexOf('K');
			while(pos >= 0){
				Peptide copy = new Peptide(p);
				copy.setPos(new int[]{pos+1});
				copy.setLinkedPos(pos+1);
//				TheoreticalSpectrum th = new TheoreticalSpectrum(copy, linkedquery.charge);
				TheoreticalSpectrum th = new LazyEvaluateLinkedSpectrum(copy, linkedquery.charge);
				//TheoreticalSpectrum th = new TheoreticalSpectrum(p);
				specList.add(th);
				pos = pep.indexOf('K', pos+1);
			}
			
			pos = pep.indexOf('G');
			while(pos >= 0){
				Peptide copy = new Peptide(p);
				copy.setPos(new int[]{pos+1});
				copy.setLinkedPos(pos+1);
//				TheoreticalSpectrum th = new TheoreticalSpectrum(copy, linkedquery.charge);
				TheoreticalSpectrum th = new LazyEvaluateLinkedSpectrum(copy, linkedquery.charge);
				//TheoreticalSpectrum th = new TheoreticalSpectrum(p);
				specList.add(th);
				pos = pep.indexOf('G', pos+1);
			}
		}
		return specList;
	}
	
	public static void testGetIonStat(){
		String file = "..\\mixture_linked\\spectrums_raw.mgf";
		SpectrumLib lib = new SpectrumLib(file, "MGF");
		LinkedPeakScoreLearner learner = new LinkedPeakScoreLearner(lib);
		learner.getIonsCount();
		//learner.getNoisesStatistics();
	}
	
	public static void testScoreFilter(){
		String libfile = "..\\mixture_linked\\linked_peptide_spectra2.mgf";
		SpectrumLib lib = new SpectrumLib(libfile, "MGF");
		//LinkedPeakScoreLearner learner = new LinkedPeakScoreLearner(lib);
		//learner.getIonsCount();
		String libfile2 = "..\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib2 = new SpectrumLib(libfile2, "MSP");
		//String mixlibfile = "..\\mixture_linked\\mixture100000.mgf";
		lib2.removeModSpectra();
		//lib2.windowFilterPeaks(6, 25);
		lib2.computeRank();
		//LPeakScoreLearner learner = new LPeakScoreLearner(lib2);
		LPeakRankBaseScorer learner = new LPeakRankBaseScorer(lib2);
		lib2 = null;
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusLinked_peptides_plusDecoy.txt";
		//String file = "..\\mixture_linked\\tempLinkedPeptides.txt";
		List<Spectrum> specList = lib.getSpectrumList();
		CandidateSpectrumLibFactory factory =  
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		LookUpSpectrumLibX lookup = new LookUpSpectrumLibX(factory.peptides, 0.5);
		factory = null;
		System.out.println("Done indexing peptides");
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		scorer.matchTolerance = 0.3;
		scorer.includeNoise = false;
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			s.windowFilterPeaks(10, 25);
			s.computePeakRank();
			List<Peak> queryPeaks = s.getTopPeaks(20);
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<String> candidates = lookup.getSpectrumByPeaks(queryPeaks, 3, s);
			String[] peptides = s.peptide.split(" & ");
			//int passedFilter = lookup.checkPassFilter(peptides[0], peptides[1].substring(0, peptides[1].length()-2), candidates);
			//int passedFilter = lookup.checkPassFilter(peptides[0], peptides[1], candidates);
			System.out.println("Query: "  + s.spectrumName + " After filter one we have: " + candidates.size() + " candidates ");
			//System.out.println("After filter correct peptide is retained?: " + passedFilter);
//			if(candidates.size() > 20000){
//				System.out.println("skipping " + s.spectrumName);
//				continue;
//			}
			List<Spectrum> candidateSpectrum = generateSpectra(candidates, s);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidateSpectrum, scorer);
			//int[] ranks = searcher.linkedRanks(s);
			//System.out.println("target peptides ranks " + ranks[0] + "\t" + ranks[1]);
			Spectrum[] topSpectra = searcher.topSpectra(s, 20);
			List<Spectrum> candidatePairs = new ArrayList<Spectrum>();
			for(int j = 0; j < topSpectra.length && topSpectra[j] != null; j++){
				//List<Peak> queryPeaks2 = s.getProjectedTopPeaks((TheoreticalSpectrum)topSpectra[j], 20);
				String p = topSpectra[j].peptide.split("\\.")[0];
				double candParentMass = LookUpSpectrumLibX.getLinkedPartnerParentmass(p, s, CrossLinker.DSP);
				//List<String> candidates2 = lookup.getSpectrum(queryPeaks, 0, candParentMass, 0.5);
				List<String> candidates2 = lookup.getSpectrumByMass(candParentMass);
				for(int k = 0; k < candidates2.size(); k++){
					String pep2 = candidates2.get(k);
					//if(!pep2.equals(peptides[1]) && !pep2.equals(peptides[0])){
					//	continue;
					//}
					Peptide p2 = new Peptide(candidates2.get(k) + ".2");
					p2.setPtmmasses(new double[]{LookUpSpectrumLibX.getLinkedOffSet(candidates2.get(k), s)});
					int pos = pep2.indexOf('K');
					while(pos >= 0){
						Peptide copy = new Peptide(p2);
						copy.setPos(new int[]{pos+1});
						Peptide o = ((TheoreticalSpectrum)topSpectra[j]).p;
						o.setCharge((short)2); //artifact from making linked peptide need to fix
						TheoreticalSpectrum th = new LazyEvaluateLinkedSpectrum(o, copy, (short)s.charge, CrossLinker.DSS.getLinkerMassOffSet());
						//TheoreticalSpectrum th = new TheoreticalSpectrum(p);
						candidatePairs.add(th);
						pos = pep2.indexOf('K', pos+1);
					}
				}
				candidates2.add(p);
				Set<String> uniques = new HashSet<String>();
//				System.out.println("filtering mass is: " + candParentMass + " looking for mass: " 
//						+ PeptideMassAnalysis.computeMolecularMass(peptides[0])
//						+ " or " + PeptideMassAnalysis.computeMolecularMass(peptides[1].substring(0, peptides[1].length()-2)));
				//System.out.println("filtering mass is: " + candParentMass + " looking for mass: " 
				//		+ PeptideMassAnalysis.computeMolecularMass(peptides[0])
				//		+ " or " + PeptideMassAnalysis.computeMolecularMass(peptides[1]));
//				passedFilter = lookup.checkPassFilter(peptides[0], peptides[1].substring(0, peptides[1].length()-2), candidates2);
//				passedFilter = lookup.checkPassFilter(peptides[0], peptides[1], candidates2);
				System.out.print("Query: "  + s.spectrumName + " With candidate " + p  + " After filter two we have: " + + candidates2.size() + " candidates ");
//				System.out.println("After filter correct peptide is retained?: " + passedFilter);
			}
			searcher = new SpectrumLibSearcher(candidatePairs, scorer);
			//searcher.topSpectrum(s);
			searcher.topLinkedSpectra(s, 1);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		
	}
	
	public static void generateMixSpectra(){
		//String libfile2 = ".\\MSPLib\\Lib\\yeast.msp";
		//SpectrumLib lib2 = new SpectrumLib(libfile2, "MSP");
		String libfile2 = "..\\mixture_linked\\human_heck_training.mgf";
		SpectrumLib lib2 = new SpectrumLib(libfile2, "MGF");
		lib2.removeModSpectra();
		lib2.createMix("..\\mixture_linked\\mixtures.mgf", 100, 0.1, 0.00001, 1.0, 3, false, 5);
//		for(int i = 0; i < 1; i++){
//			SpectrumLib mixture = lib2.createRandomMix(5000, 0.1, 0.001, 1.0, 3, false);
//			System.out.println("Generated library size: " + mixture.getAllSpectrums().size());
//			mixture.printLibToFile("..\\mixture_linked\\mixture.mgf_part"+(i+1), mixture);
//		}
	}
	
	public static void main(String[] args){
		//testGetIonStat();
		//testScoreFilter();
		generateMixSpectra();
	}

}
