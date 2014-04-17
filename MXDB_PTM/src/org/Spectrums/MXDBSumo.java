package org.Spectrums;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import mixdb.ArraySpectrumComparator;
import mixdb.LinkedScorerAdapter;
import mixdb.LinkedTheoSpectrumFactory;
import SeqDB.DatabaseIndexer;
import SeqDB.PeptideUtils;

/**
 * Perform MXDBSearch for SUMOylated peptides
 * @author Jian Wang
 *
 */
public class MXDBSumo {
	public double windowWidth = 25;
	public int topPeaksKept = 10;
	//candidate filtering
	public int mode = 0;
	public int fragmentMode = 0;
	public double parentMassTolerance = 1;
	public double fragmentMassTolerance = 0.5;
	//input and training file
	public String queryFile;
	public String training = "..\\MSPLib\\Lib\\ecoli.msp";
	public String mixtureTraining = "";
	public int minScan = 0;
	public int maxScan = 0;
	public String outFile = "";
	public String[] ptms = new String[]{"Z", "QQQTGG", "Q-17.027QQTGG"};
	public String[] Ncuts=new String[]{"K", "R"};
	public String[] Ccuts=new String[]{"K", "R"};
	public int miscleaves=1;
	
	public MXDBSumo(){
		
	}
		
	public void search(String queryFile, String peptideFile){
		Iterator<Spectrum> iter= null;
		if(queryFile.endsWith(".mgf")){
			iter = new LargeSpectrumLibIterator(queryFile);
		}
		if(queryFile.endsWith(".mzXML")){
			iter = new MZXMLReader(queryFile);
		}
    	LookUpSpectrumLibXX lookup = new LookUpSpectrumLibXX();
    	//peptideFile = "../mixture_linked/database/lib_sumo_peptides_plusHuman_plusDecoy.fasta";
    	DatabaseIndexer db = new DatabaseIndexer(peptideFile, this.Ncuts, this.Ccuts, this.miscleaves, 0.01);
    	//DatabaseIndexer db = new DatabaseIndexer(peptideFile, 0.01);
    	this.queryFile = queryFile;
    	BufferedWriter bw = Utils.FileIOUtils.initOutputStream(this.outFile);
    	lookup.setParentMassTolerance(this.parentMassTolerance);
    	lookup.setMinMatchedPeak(-1);
    	lookup.setMinContinuousMatch(-1);
    	lookup.setMinCharge(1);
    	lookup.setMaxCharge(1);
    	lookup.loadPeptidesFromFileLite(peptideFile);
    	lookup.setToleranceMode(this.mode);
    	System.out.println("start searching");
    	long start = (new GregorianCalendar()).getTimeInMillis();
    	SimpleProbabilisticScorer scorer1 = (SimpleProbabilisticScorer)SpectrumUtil.getLPeakRankBaseScorer(training);
    	scorer1.matchToleranceMode = this.fragmentMode;
    	scorer1.matchTolerance = this.fragmentMassTolerance;
    	//SimpleProbabilisticScorer filter = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedPeptideSingleScorer(mixtureTraining);
    	//SimpleProbabilisticScorer scorer2 = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedSUMOScorer(mixtureTraining);
    	SimpleProbabilisticScorer scorer2 = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedSUMOScorer(mixtureTraining);
    	LinkedPeptidePeakScoreLearner peakscorer = LinkedPeptidePeakScoreLearner.loadComparator(mixtureTraining);
    	scorer2.matchToleranceMode = this.fragmentMode;
    	scorer2.matchTolerance = this.fragmentMassTolerance;
    	//SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(mixtureTraining);
//   	 	((LinkedPeptidePeakScoreLearner)((MixtureSpectrumScorer)scorer2).getComp()).generalLinkedMode=true;
		LinkedScorerAdapter adaptor = new LinkedScorerAdapter(LinkedTheoSpectrumFactory.linkedTypeMap, 
				LinkedTheoSpectrumFactory.linkedIonMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		arryScorer.setScoreTable(adaptor.getTable());
		arryScorer.setErrorTable(adaptor.getErrorsTable());
    	
		List<Peak> pList = new ArrayList();
		//creat sumo
		List<Peptide> peptides2 = new ArrayList();
		System.out.println("Number of PTM tag consider: " + this.ptms.length);
		for(int i =0; i < this.ptms.length; i++){
			peptides2.add(new Peptide(ptms[i], 1));
		}
		Mass.DSSLINKER_MASS = -1*Mass.WATER;
		CrossLinker SUMO = new CrossLinker(Mass.WATER*-1, new int[]{CrossLinker.ANYPOSITION}, new int[]{CrossLinker.CTERM, CrossLinker.ANYPOSITION}, new char[]{'K'}, new char[]{'G','Z'});
    	//TheoreticalSpectrum.deconvolutedMode = true;
		for(;iter.hasNext();){
    		Spectrum s = iter.next();
    		//System.out.println("Scan: " + s.scanNumber);
    		if(s.scanNumber < this.minScan || s.scanNumber > this.maxScan 
    					|| s.charge >= 6 || s.charge < 2){
    			continue;
    		}
    		//s.charge = (s.charge+1)*-1;
    		s.removePrecursors(0.5);
    		//s = new DeconvolutedSpectrum(s);
    		s.windowFilterPeaks(this.topPeaksKept, this.windowWidth);
    		//SpectrumUtil.removeSUMOPeaks(s);
    		//s.windowFilterPeaks(8, 25);
    		s.computePeakRank();
    		System.out.println("Searching " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge +"\t" + s.getPeak().size());
			List<Peptide> linkedpeps = new ArrayList();
    		int minCharge = s.charge > 0 ? s.charge : 2;
    		int maxCharge = s.charge > 0 ? s.charge : 4;
			for(int charge = minCharge; charge <=maxCharge; charge++){
    			for(int i = 0; i < peptides2.size(); i++){
    				Peptide p2 = peptides2.get(i);
    				double lookupmass = s.parentMass*charge - Mass.PROTON_MASS*(charge-1) - p2.getParentmass() - SUMO.getLinkerMassOffSet() + Mass.PROTON_MASS;
    				//System.out.println("lookup mass " + lookupmass);
    				//List<PeptideLite> candidates = lookup.getCandidatePeptide(lookupmass, charge, pList);
    				//List<Peptide> peps = LinkedPeakScoreLearner.generatePeptides(candidates);
    				lookupmass = lookupmass - Mass.WATER - Mass.PROTON_MASS;
    				//List<PeptideLite> candidates = db.getPeptidesWithC13(lookupmass - (lookupmass*this.parentMassTolerance)/1000000, lookupmass + lookupmass*this.parentMassTolerance/1000000, 0.05,1);
    				//List<Peptide> peps = PeptideUtils.generatePeptide(candidates, db.getSeq(), 0.05);
    				List<Peptide> peps = db.getPeptidesFull(lookupmass - ((lookupmass*this.parentMassTolerance)/1000000.0), 
    						lookupmass + (lookupmass*this.parentMassTolerance)/1000000, 1);
    				lookupmass = lookupmass + Mass.WATER + Mass.PROTON_MASS;
    				
    				double massoffset = lookupmass - 128.0949 + SUMO.getLinkerMassOffSet() - Mass.PROTON_MASS; //lysine mass
    				if(massoffset > 0){
    					peps.add(new Peptide("Z"+massoffset+"K.1"));
    				}
    				double lookupmass1 = s.parentMass*charge - Mass.PROTON_MASS*(charge-1) - p2.getParentmass() - SUMO.getLinkerMassOffSet() + Mass.PROTON_MASS - 42.010565;
    				//List<PeptideLite> candidates1 = lookup.getCandidatePeptide(lookupmass1, charge, pList);
    				//List<Peptide> peps1 = LinkedPeakScoreLearner.generatePeptides(candidates1);
    				
    				lookupmass1 = lookupmass1 - Mass.WATER - Mass.PROTON_MASS;
    				//List<PeptideLite> candidates1 = db.getPeptidesWithC13(lookupmass1 - (lookupmass1*this.parentMassTolerance)/1000000, lookupmass1 + lookupmass1*this.parentMassTolerance/1000000, 0.05,1);
    				//List<Peptide> peps1 = PeptideUtils.generatePeptide(candidates1, db.getSeq(), 0.05);
    				List<Peptide> peps1 = db.getPeptidesFull(lookupmass1 - (lookupmass1*this.parentMassTolerance)/1000000, 
    						lookupmass1 + (lookupmass1*this.parentMassTolerance)/1000000, 1);
    				
    				peps1 = Peptide.insertPTM(peps1, 42.010565, 1, 1);
    				peps.addAll(peps1);
    				double lookupmass2 = s.parentMass*charge - Mass.PROTON_MASS*(charge-1) - p2.getParentmass() - SUMO.getLinkerMassOffSet() + Mass.PROTON_MASS - 15.994915;
    				//List<PeptideLite> candidates2 = lookup.getCandidatePeptide(lookupmass2, charge, pList);
    				//List<Peptide> peps2 = LinkedPeakScoreLearner.generatePeptides(candidates2);
    				
    				lookupmass2 = lookupmass2 - Mass.WATER - Mass.PROTON_MASS;
    				//List<PeptideLite> candidates2 = db.getPeptidesWithC13(lookupmass2 - (lookupmass2*this.parentMassTolerance)/1000000, lookupmass2 + (lookupmass2*this.parentMassTolerance/1000000), 0.05,1);
    				//List<Peptide> peps2 = PeptideUtils.generatePeptide(candidates2, db.getSeq(), 0.05);
    				List<Peptide> peps2 = db.getPeptidesFull(lookupmass2 - (lookupmass2*this.parentMassTolerance)/1000000, 
    						lookupmass2 + (lookupmass2*this.parentMassTolerance)/1000000, 1);
    				
    				peps2 = Peptide.insertPTM(peps2, 15.994915, new char[]{'M'}, 1);
    				//peps.addAll(peps2);
    				//double lookupmass3 = s.parentMass*charge - Mass.PROTON_MASS*(charge-1) - p2.getParentmass() - SUMO.getLinkerMassOffSet() + Mass.PROTON_MASS - 8.041;
    				//List<PeptideLite> candidates3 = lookup.getCandidatePeptide(lookupmass3, charge, pList);
    				//List<Peptide> peps3 = LinkedPeakScoreLearner.generatePeptides(candidates3);
    				//peps3 = Peptide.insertPTM(peps3, 8.041, new char[]{'K'}, 1);
    				//peps.addAll(peps3);
    				linkedpeps.addAll(CandidateSpectrumLibFactory.getCrossLinkPeptides(peps, p2, SUMO, charge, charge));
    			}
    		}
    		List<Spectrum> sList = new ArrayList();
    		for(Iterator it = linkedpeps.iterator(); it.hasNext();){
    			LinkedPeptide cand = (LinkedPeptide)it.next();
    			if(cand.peptides[0].getLinkedPos() == cand.peptides[0].getPeptide().length()){ //no cterm sumo linkage
    				continue;
    			}
    			LazyEvaluateLinkedSpectrum t = new LazyEvaluateLinkedSpectrum(((LinkedPeptide)cand).peptides[0], 
						((LinkedPeptide)cand).peptides[1], cand.getCharge(), SUMO.getLinkerMassOffSet());
    			//System.out.println("candidate is: " + cand);
    			sList.add(t);
    		}
    		System.out.println("Scoring " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge + "\tnumber of candidates: " + sList.size());
    		//SpectrumLibSearcher searcher = new SpectrumLibSearcher(sList, scorer2);
    		//SpectrumLibSearcher searcher = new SpectrumLibSearcher(sList, scorer2);
    		SpectrumLibSearcher searcher = new SpectrumLibSearcher(sList, arryScorer);
			searcher.setSingleScorer(scorer1);
			if(this.fragmentMode==0){
				searcher.matchTolerance = this.fragmentMassTolerance;
			}else{
				searcher.matchTolerance = this.fragmentMassTolerance/1000;
			}
			//searcher.topLinkedSpectra(s, 10);
			searcher.queryFile = this.queryFile;
			searcher.bw = bw;
			searcher.topArrayCandidates(s, 1);
			try{
				bw.flush();
			}catch(IOException ioe){
				
			}
    	}
		Utils.FileIOUtils.finishOutput(bw);
 	}
	
	public static void testMXDBSUMO(String inFile){
		Map<String, String> arguments= Utils.FileIOUtils.createTableFromFile(inFile, 0, 1);
		MXDBSumo mxdb = new MXDBSumo();
		String[] prefixes = arguments.get("PrefixIons").split(",");
		String[] suffixes = arguments.get("SuffixIons").split(",");
		TheoreticalSpectrum.prefixIons = prefixes;
		TheoreticalSpectrum.suffixIons = suffixes;
		String[] ptms = arguments.get("PTMTag").split(",");
		mxdb.ptms = ptms;
		mxdb.Ncuts = arguments.get("NtermCut").split(",");
		mxdb.Ccuts = arguments.get("CtermCut").split(",");
		mxdb.miscleaves = Integer.parseInt(arguments.get("Miscleaves"));
		mxdb.mixtureTraining = arguments.get("MixtureModel");
		mxdb.training = arguments.get("SingleModel");
		mxdb.parentMassTolerance = Double.parseDouble(arguments.get("ParentMassTolerance"));
		mxdb.fragmentMassTolerance = Double.parseDouble(arguments.get("FragmentMassTolerance"));
		mxdb.topPeaksKept = Integer.parseInt(arguments.get("TopPeaksKept"));
		mxdb.windowWidth = Double.parseDouble(arguments.get("WindowWidth"));
		String queryFile = arguments.get("QueryFile");
		String peptideFile = arguments.get("PeptideFile");
		mxdb.mode = Integer.parseInt(arguments.get("ToleranceMode"));
		mxdb.fragmentMode = Integer.parseInt(arguments.get("FragmentMode"));
		mxdb.minScan = Integer.parseInt(arguments.get("MinScan"));
		mxdb.maxScan = Integer.parseInt(arguments.get("MaxScan"));
		mxdb.outFile = arguments.get("OutputFile");
		mxdb.search(queryFile, peptideFile);
	}
	
	public static void main(String[] args){
		args[0] = "..\\mixture_linked\\MXDBSumo_inputs.txt";
		testMXDBSUMO(args[0]);
	}
}
