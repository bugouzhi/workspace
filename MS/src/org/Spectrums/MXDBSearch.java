package org.Spectrums;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
	public char linkerSite = 'K';
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
    	LookUpSpectrumLibXX lookup = new LookUpSpectrumLibXX();
    	lookup.setMinCharge(1);
    	lookup.setMaxCharge(1);
    	//List<String> peptides = Utils.FileIOUtils.createListFromFile(peptideFile);
    	//LookUpSpectrumLibX pLookup = new LookUpSpectrumLibX(peptides, this.parentMassTolerance, LookUpSpectrumLibX.PARENTMODE);
    	lookup.setMinMatchedPeak(this.minMatchedPeak);
    	lookup.setMinContinuousMatch(this.minContinuousMatch);
    	lookup.loadPeptidesFromFile(peptideFile);
    	System.out.println("start searching");
    	long start = (new GregorianCalendar()).getTimeInMillis();
    	SimpleProbabilisticScorer scorer1 = (SimpleProbabilisticScorer)SpectrumUtil.getLPeakRankBaseScorer(training);
    	SimpleProbabilisticScorer filter = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedPeptideSingleScorer(mixtureTraining);   	
//    	SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(mixtureTraining);
    	SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(mixtureTraining);
    	//((LinkedPeptidePeakScoreLearner)((MixtureSpectrumScorer)scorer2).getComp()).generalLinkedMode=true;
//    	List<Spectrum> specList = reader.readAllMS2Spectra();
//    	Iterator<Spectrum> iter = specList.iterator();
    	Mass.DSSLINKER_MASS = this.linkerMass;
    	for(;iter.hasNext();){
    		Spectrum s = iter.next();
    		//s.charge = 4;
    		if(s.scanNumber < this.minScan || s.scanNumber > this.maxScan || s.charge >= 8 || s.charge < 2){
    			continue;
    		}
    		//s.charge = (s.charge+1)*-1;
    		//s.charge = 3;
    		s.windowFilterPeaks(this.topPeaksKept, this.windowWidth);
    		s.removePrecursors(0.5);
    		s.computePeakRank();
    		List<Peak> pList = s.getTopPeaks(this.topFilteringPeak);
    		System.out.println("Searching " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
    		List<PeptideLite> candidates = lookup.getCandidatePeptide(s.parentMass, s.charge, pList);
    		String[] targetPeptides = s.peptide.split("--");
			//int passedFilter = lookup.checkPassFilter(targetPeptides[0], targetPeptides[1], candidates);
			//System.out.println("Query: "  + s.spectrumName + " after filter one we have: " + candidates.size() + " candidates\tAfter filter correct peptide is retained?: " + passedFilter);
			if(true){
				//continue;
			}
    		List<Peptide> peps = LinkedPeakScoreLearner.generatePeptides(candidates);
    		List<Peptide> candPepWithMod = new ArrayList();
    		candPepWithMod.addAll(peps);
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 156.0680, new char[]{'K'}, 1));
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 57, 1, 1));
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 15.995, new char[]{'M'}, 1));
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 42.010565, 1, 1));
    		List<Peptide> linkedPeps = LinkedPeakScoreLearner.generateLinkedPeptides(candPepWithMod, s, this.linkerSite);
    		System.out.println(s.spectrumName + " has candidates: " + linkedPeps.size());
    		List<Spectrum> candidateSpectrum = LinkedPeakScoreLearner.generateSpectra3(linkedPeps, s);
   			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidateSpectrum, filter);
			if(targetPeptides.length == 2){
				int[] ranks = searcher.linkedRanks(s);
				System.out.println(s.spectrumName + " target peptides ranks " + ranks[0] + "\t" + ranks[1]);
			}
			if(targetPeptides.length > 0){
				//continue;
			}
			Spectrum[] topSpectra = searcher.topSpectra(s, this.topFirstPassMatch);

			List<Spectrum> candidatePairs = new ArrayList<Spectrum>();
			for(int j = 0; j < topSpectra.length && topSpectra[j] != null; j++){
				String p = topSpectra[j].peptide.split("\\.")[0];
				Peptide p1 = ((TheoreticalSpectrum)topSpectra[j]).p;
				Peptide o = new Peptide(p1);
				char aaMatch = 'B';
				if(p.charAt(o.getLinkedPos()-1) == 'K'){
					aaMatch = 'K';
				}
				if(p.charAt(o.getLinkedPos()-1) == 'G'){
					aaMatch = 'G';
				}
				aaMatch = this.linkerSite;
				double candParentMass = LookUpSpectrumLibX.getLinkedPartnerParentmass(p1, s, CrossLinker.DSS);
				//System.out.println("pep1: " + p1 + "\tlookup: " + candParentMass);
				
				//%%%List<String> candidates2_old = pLookup.getSpectrumByMass(candParentMass);
				//%%%candidates2.addAll(pLookup.getSpectrumByMass(candParentMass-1.0));
				
				lookup.setParentMassTolerance(this.parentMassTolerance);
				List<PeptideLite> candidates2lite = lookup.getCandidatePeptide(candParentMass+Mass.PROTON_MASS, 1, pList);
				//System.out.println("peptide2_old size: " + candidates2_old.size());
				//System.out.println("peptide2 size: " + candidates2lite.size());
	    		List<Peptide> candidates2 = LinkedPeakScoreLearner.generatePeptides(candidates2lite);
	    		//System.out.println("peptide2 size: " + candidates2.size());
	    		
				for(int k = 0; k < candidates2.size(); k++){
					//System.out.println("spect: " + s.parentMass + "\t" + s.charge);
					//System.out.println("pairing1: " + o + "\t" + o.getParentmass()+ "\twith mass: " + candParentMass);
					//%%%%String pep2 = candidates2.get(k);
					//%%%%Peptide p2 = new Peptide(candidates2.get(k) + ".2");
					Peptide p2 = new Peptide(candidates2.get(k));
					String pep2 = p2.getPeptide();
					double offset2 = LookUpSpectrumLibX.getLinkedOffSet(candidates2.get(k), s);
					//p2.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(o.getPeptide())-Mass.WATER});
					//o.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(candidates2.get(k))-Mass.WATER});
					int pos = pep2.indexOf(aaMatch);
					//System.out.println("offset2: " + offset2);
					while(pos >= 0 && pos != pep2.length()-1){
						Peptide copy = new Peptide(p2);
						copy.insertPTM(pos+1, offset2);
						copy.setLinkedPos(pos+1);
						//o.setCharge((short)2); //artifact from making linked peptide need to fix
						//System.out.println("peptide1: " + o + "\t"+ o.getCharge() + "\tPeptide2: "  + copy +"\t" + copy.getCharge());
						//TheoreticalSpectrum th = new TheoreticalSpectrum(o, copy, (short)s.charge, true);
						LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(o, copy, (short)s.charge);
						//System.out.println("peptide generated is: " + th.peptide);
						candidatePairs.add(th);
						pos = pep2.indexOf(aaMatch, pos+1);
					}
				}
				
				double modMass = 42.010565;//15.995;
				char modAA = 'M';
				double candParentMassMod = candParentMass - modMass;
				//%%%List<String> candidates2Mod = pLookup.getSpectrumByMass(candParentMassMod);
				//%%%candidates2Mod.addAll(pLookup.getSpectrumByMass(candParentMassMod-1.0));
				
				List<PeptideLite> candidates2liteMod = lookup.getCandidatePeptide(candParentMassMod+Mass.PROTON_MASS, 1, pList);
	    		List<Peptide> candidates2Mod = LinkedPeakScoreLearner.generatePeptides(candidates2liteMod);
	    		
				for(int k = 0; k < candidates2Mod.size(); k++){
					//System.out.println("pairing: " + o + "\twith mass: " + candParentMassMod);
					
					//%%%String pep2 = candidates2Mod.get(k);
					//%%%Peptide p2 = new Peptide(candidates2Mod.get(k) + ".2");	
				
					Peptide p2 = new Peptide(candidates2Mod.get(k));
					String pep2 = p2.getPeptide();
					int modInd = p2.getPeptide().indexOf('M');
					modInd = 0;
					if(modInd < 0){
						continue;
					}
					p2.insertPTM(modInd+1, modMass);
					double offset2 = LookUpSpectrumLibX.getLinkedOffSet(candidates2Mod.get(k), s) - modMass;
					//System.out.println("offset2: " + offset2);
					//p2.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(o.getPeptide())-Mass.WATER});
					//o.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(candidates2.get(k))-Mass.WATER});
					int pos = pep2.indexOf(aaMatch);
					while(pos >= 0 && pos != pep2.length()-1 && pos!= 0){  //we disallow double PTM at the same residue position
						Peptide copy = new Peptide(p2);
						int posM = pep2.indexOf(modAA);
						//copy.insertPTM(posM+1, modMass);
						copy.insertPTM(pos+1, offset2);
						copy.setLinkedPos(pos+1);
						//o.setCharge((short)2); //artifact from making linked peptide need to fix
						//System.out.println("peptide1: " + o + "\t"+ o.getCharge() + "\tPeptide2: "  + copy +"\t" + copy.getCharge());
						//TheoreticalSpectrum th = new TheoreticalSpectrum(o, copy, (short)s.charge, true);
						LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(o, copy, (short)s.charge);
						//System.out.println("peptide generated is: " + th.peptide);
						candidatePairs.add(th);
						pos = pep2.indexOf(aaMatch, pos+1);
					}
				}
				//added a dummy cross-link where the 2nd peptide is unspecified
				o = new Peptide(p1);
				Peptide p2 = new Peptide("ZZ.2");
				double offset1 = p1.getPtmmasses()[p1.getPtmmasses().length-1] - Mass.DSSLINKER_MASS - Mass.WATER; //need to modify, now assume linked peptide is alway last mod
				p2.insertPTM(1, offset1);
				double offset2 = LookUpSpectrumLibX.getLinkedOffSet(p2, s);
				//System.out.println("Peptide is: " + p1 + "\toffset:\t" + offset1 + "\t" + offset2);
				p2.insertPTM(2, offset2);
				p2.setLinkedPos(2);
				LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(o, p2, (short)s.charge);
				//System.out.println("peptide generated isa: " + th.peptide);
				//candidatePairs.add(th);
				//candidates2.add(p);
				Set<String> uniques = new HashSet<String>();
			}
			System.out.println(s.spectrumName + " has candidates pairs: " + candidatePairs.size());
			lookup.setParentMassTolerance(3000.0);
			searcher = new SpectrumLibSearcher(candidatePairs, scorer2);
			searcher.matchTolerance = this.fragmentMassTolerance;
			searcher.setSingleScorer(scorer1);
			searcher.topSpectrum(s);
			searcher.topLinkedSpectra(s, 1, -100);
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
		mxdb.mixtureTraining = arguments.get("MixtureTraining");
		mxdb.training = arguments.get("Training");
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
		mxdb.linkerMass = Double.parseDouble(arguments.get("LinkerMass"));
		mxdb.linkerSite = arguments.get("LinkerSite").charAt(0);
		mxdb.search(queryFile, peptideFile);
	}
	
	public static void main(String[] args){
		args[0] = "..\\mixture_linked\\MXDB_inputs2.txt";
		testMXDBSearch(args[0]);
	}
}
