package org.Spectrums;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import SeqDB.DatabaseIndexer;

import mixdb.ArraySpectrumComparator;
import mixdb.LinkedScorerAdapter;
import mixdb.LinkedTheoSpectrumFactory;

/**
 * Search for mixture and crosslinked peptide
 * @author bugouzhi
 *
 */
public class MXDBSearch {
	//query filtering
	public double windowWidth = 25;
	public int topPeaksKept = 10;
	//candidate filtering
	public double parentMassTolerance = 1;
	public int minMatchedPeak = 3; //for filtering
	public int minContinuousMatch = 0; // for filtering
	public int topFilteringPeak = 20; //for filtering
	//public int minContinuousMatchedPeak = 1;
	//searching
	public int topFirstPassMatch = 25; //top match to retain in the first pass
	public double fragmentMassTolerance = 0.5;
	//input and training file
	public String queryFile;
	public String training = "..\\MSPLib\\Lib\\ecoli.msp";
	public String mixtureTraining = "";
	public int minScan = 0;
	public int maxScan = 0;
	public double linkerMass = 138.0680;
	public double[] linkerMasses;
	public char linkerSite = 'K';
	public String[] Ncuts=new String[]{"K", "R"};
	public String[] Ccuts=new String[]{"K", "R"};
	public int miscleaves=1;
	public int numC13 = 1;
	public double CysProt = 57.021410;
	public String outFile="";
	
	
	public MXDBSearch(){
		
	}
		
	public void search(String queryFile, String peptideFile){
		Iterator<Spectrum> iter= null;
		if(queryFile.endsWith(".mgf")){
			iter = new LargeSpectrumLibIterator(queryFile);
		}
		if(queryFile.endsWith(".mzXML")){
			iter = new MZXMLReader(queryFile);
		}
		this.queryFile = queryFile;
    	System.out.println("start searching");
    	long start = (new GregorianCalendar()).getTimeInMillis();
    	SimpleProbabilisticScorer scorer1 = (SimpleProbabilisticScorer)SpectrumUtil.getLPeakRankBaseScorer(training);
    	SimpleProbabilisticScorer filter = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedPeptideSingleScorer(mixtureTraining);   	
//    	SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(mixtureTraining);
    	SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(mixtureTraining);
    	//((LinkedPeptidePeakScoreLearner)((MixtureSpectrumScorer)scorer2).getComp()).generalLinkedMode=true;
      	LinkedPeptidePeakScoreLearner peakscorer = LinkedPeptidePeakScoreLearner.loadComparator(mixtureTraining);      	 
    	LinkedScorerAdapter adaptor = new LinkedScorerAdapter(LinkedTheoSpectrumFactory.linkedTypeMap, 
				LinkedTheoSpectrumFactory.linkedIonMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		arryScorer.setScoreTable(adaptor.getTable());
		arryScorer.setErrorTable(adaptor.getErrorsTable());
		//setting mass parameters
		Mass.addFixMod('C', this.CysProt-57.021410);
		Mass.DSSLINKER_MASS = this.linkerMass;
    
		
		// We use to indexer here because the first stage, which is essentially a blind search, has very wide tolerance
		//it will be very inefficient to use a index table with very small resolution
		DatabaseIndexer db = new DatabaseIndexer(peptideFile, this.Ncuts, this.Ccuts, this.miscleaves, 100);
		DatabaseIndexer db2 = new DatabaseIndexer(peptideFile, this.Ncuts, this.Ccuts, this.miscleaves, 0.01);
		BufferedWriter bw = Utils.FileIOUtils.initOutputStream(this.outFile);
		//System.out.println("this linker " + this.linkerMass);
    	CrossLinker linker = new CrossLinker(this.linkerMasses, new char[]{'K'});
    	for(;iter.hasNext();){
    		Spectrum s = iter.next();
    		//s.charge = 4;
    		if(s.scanNumber < this.minScan || s.scanNumber > this.maxScan || s.charge >= 8 || s.charge < 3){
    			continue;
    		}
    		//s.charge = (s.charge+1)*-1;
    		s.windowFilterPeaks(this.topPeaksKept, this.windowWidth);
    		s.removePrecursors(0.5);
    		s.computePeakRank();
    		//System.out.println("Searching " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
    		String[] targetPeptides = s.peptide.split("--");
			//int passedFilter = lookup.checkPassFilter(targetPeptides[0], targetPeptides[1], candidates);
			//System.out.println("Query: "  + s.spectrumName + " after filter one we have: " + candidates.size() + " candidates\tAfter filter correct peptide is retained?: " + passedFilter);
			
    		//List<Peptide> peps = LinkedPeakScoreLearner.generatePeptides(candidates);
    		
    		List<Peptide> peps = db.getPeptidesFull(500, s.parentMass*s.charge - s.charge*Mass.PROTON_MASS -Mass.WATER -200,
					1);
    		
    		List<Peptide> candPepWithMod = new ArrayList();
    		candPepWithMod.addAll(peps);
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 156.0680, new char[]{'K'}, 1));
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 57, 1, 1));
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 15.995, new char[]{'M'}, 1));
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 42.010565, 1, 1));
    		List<Peptide> linkedPeps =LinkedPeptide.generateLinkedPeptides(candPepWithMod, s, this.linkerSite);
    		System.out.println(s.spectrumName + " has candidates: " + linkedPeps.size());
    		List<Spectrum> candidateSpectrum = LinkedPeakScoreLearner.generateSpectra3(linkedPeps, s);
   			
    		//SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidateSpectrum, filter);
    		SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidateSpectrum, arryScorer);
    		if(targetPeptides.length == 2){
				//int[] ranks = searcher.linkedRanks(s);
				//System.out.println(s.spectrumName + " target peptides ranks " + ranks[0] + "\t" + ranks[1]);
			}
    		
			Spectrum[] topSpectra = searcher.topArrayCandidates(s, this.topFirstPassMatch, false);
			if(true){
				//continue;
			}
			List<Spectrum> candidatePairs = new ArrayList<Spectrum>();
			for(int j = 0; j < topSpectra.length && topSpectra[j] != null; j++){
				String p = topSpectra[j].peptide.split("\\.")[0];
				Peptide p1 = ((TheoreticalSpectrum)topSpectra[j]).p;
				//System.out.println("Fasta protein " + p1.getFastaseq());
				Peptide o = new Peptide(p1);
				char aaMatch = this.linkerSite;
				double[] candParentMass = LookUpSpectrumLibX.getLinkedPartnerParentmass(p1, s, linker);
				for(int k = 0; k < candParentMass.length; k++){
					List<Peptide> cands = db2.getPeptidesFull(candParentMass[k] - this.parentMassTolerance - Mass.WATER, 
	    					candParentMass[k] + this.parentMassTolerance - Mass.WATER, 
	    					this.numC13);
					//System.out.println(o + ": " + candParentMass[k] + "\t" + cands.size());
					for(int l = 0; l < cands.size(); l++){
						Peptide p2 = new Peptide(cands.get(l));
						String pep2 = p2.getPeptide();
						double offset2 = LookUpSpectrumLibX.getLinkedOffSet(cands.get(l), s);
						int pos = pep2.indexOf(aaMatch);
						//System.out.println("offset2: " + offset2);
						while(pos >= 0 && pos != pep2.length()-1){
							Peptide copy = new Peptide(p2);
							copy.insertPTM(pos+1, offset2);
							copy.setLinkedPos(pos+1);
							LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(o, copy, (short)s.charge, linker.getLinkerMassOffSets()[k]);
							candidatePairs.add(th);
							pos = pep2.indexOf(aaMatch, pos+1);
						}
					}
				}
				
				double modMass = 42.010565;//15.995;
				char modAA = 'M';
				
				double[] candParentMassMod = ArrayUtils.shift(candParentMass, -1*modMass);
				for(int k = 0; k < candParentMass.length; k++){
					List<Peptide> cands = db2.getPeptidesFull(candParentMassMod[k] - this.parentMassTolerance - Mass.WATER, 
	    					candParentMassMod[k] + this.parentMassTolerance - Mass.WATER,
	    					this.numC13);
				    		
					for(int l = 0; l < cands.size(); l++){	
						Peptide p2 = new Peptide(cands.get(l));
						String pep2 = p2.getPeptide();
						int modInd = p2.getPeptide().indexOf('M');
						modInd = 0;
						if(modInd < 0){
							continue;
						}
						p2.insertPTM(modInd+1, modMass);
						double offset2 = LookUpSpectrumLibX.getLinkedOffSet(cands.get(l), s) - modMass;
						int pos = pep2.indexOf(aaMatch);
						while(pos >= 0 && pos != pep2.length()-1 && pos!= 0){  //we disallow double PTM at the same residue position
							Peptide copy = new Peptide(p2);
							int posM = pep2.indexOf(modAA);
							//copy.insertPTM(posM+1, modMass);
							copy.insertPTM(pos+1, offset2);
							copy.setLinkedPos(pos+1);
							LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(o, copy, (short)s.charge, linker.getLinkerMassOffSets()[k]);
							candidatePairs.add(th);
							pos = pep2.indexOf(aaMatch, pos+1);
						}
					}
				}
				//added a dummy cross-link where the 2nd peptide is unspecified
				o = new Peptide(p1);
				Peptide p2 = new Peptide("ZZ.2");
				p2.setFastaseq(o.getFastaseq());
				p2.setBeginIndex(1); //just some arbitrary index
				double offset1 = p1.getPtmmasses()[p1.getPtmmasses().length-1] - Mass.DSSLINKER_MASS - Mass.WATER; //need to modify, now assume linked peptide is alway last mod
				p2.insertPTM(1, offset1);
				double offset2 = LookUpSpectrumLibX.getLinkedOffSet(p2, s);
				//System.out.println("Peptide is: " + p1 + "\toffset:\t" + offset1 + "\t" + offset2);
				p2.insertPTM(2, offset2);
				p2.setLinkedPos(2);
				//LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(o, p2, (short)s.charge, Mass.DSSLINKER_MASS);
				//System.out.println("peptide generated isa: " + th.peptide);
				//candidatePairs.add(th);
				//candidates2.add(p);
				Set<String> uniques = new HashSet<String>();
			}
			//lookup.setParentMassTolerance(3000.0);
			searcher = new SpectrumLibSearcher(candidatePairs, scorer2);
			searcher.adjustPrintOrder = true;
			System.out.println("We gathered pairs: " + candidatePairs.size());
			searcher.matchTolerance = this.fragmentMassTolerance;
			searcher.setSingleScorer(scorer1);
			//searcher.topSpectrum(s);
			searcher.queryFile = this.queryFile;
			searcher.bw = bw;
			searcher.topLinkedSpectra(s, 1);
    	}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMXDBSearch(String inFile){
		Map<String, String> arguments= Utils.FileIOUtils.createTableFromFile(inFile, 0, 1);
		MXDBSearch mxdb = new MXDBSearch();
		String[] prefixes = arguments.get("PrefixIons").split(",");
		String[] suffixes = arguments.get("SuffixIons").split(",");
		TheoreticalSpectrum.prefixIons = prefixes;
		TheoreticalSpectrum.suffixIons = suffixes;
		mxdb.mixtureTraining = arguments.get("MixtureModel");
		mxdb.training = arguments.get("SingleModel");
		mxdb.Ncuts = arguments.get("NtermCut").split(",");
		mxdb.Ccuts = arguments.get("CtermCut").split(",");
		mxdb.miscleaves = Integer.parseInt(arguments.get("Miscleaves"));
		mxdb.parentMassTolerance = Double.parseDouble(arguments.get("ParentMassTolerance"));
		mxdb.fragmentMassTolerance = Double.parseDouble(arguments.get("FragmentMassTolerance"));
		mxdb.topFirstPassMatch = Integer.parseInt(arguments.get("TopFirstPassMatch"));
		mxdb.topPeaksKept = Integer.parseInt(arguments.get("TopPeaksKept"));
		mxdb.windowWidth = Double.parseDouble(arguments.get("WindowWidth"));
		mxdb.minMatchedPeak = -1;
		mxdb.minContinuousMatch = -1;
		String queryFile = arguments.get("QueryFile");
		String peptideFile = arguments.get("PeptideFile");
		mxdb.minScan = Integer.parseInt(arguments.get("MinScan"));
		mxdb.maxScan = Integer.parseInt(arguments.get("MaxScan"));
		mxdb.linkerMasses = CommandArgument.getDoulbArray(arguments.get("LinkerMass"), ",");
		mxdb.linkerSite = arguments.get("LinkerSite").charAt(0);
		mxdb.outFile = arguments.get("OutputFile");
		mxdb.numC13 = Integer.parseInt(arguments.get("NumC13"));
		//mxdb.CysProt = Double.parseDouble(arguments.get("CysProt"));
		mxdb.search(queryFile, peptideFile);
	}
	
	public static void main(String[] args){
		//args[0] = "..\\mixture_linked\\MXDB_inputs2.txt";
		testMXDBSearch(args[0]);
	}
}
