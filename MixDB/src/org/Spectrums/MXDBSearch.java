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
    	List<String> peptides = Utils.FileIOUtils.createListFromFile(peptideFile);
    	LookUpSpectrumLibX pLookup = new LookUpSpectrumLibX(peptides, this.parentMassTolerance, LookUpSpectrumLibX.PARENTMODE);
    	lookup.setMinMatchedPeak(this.minMatchedPeak);
    	lookup.setMinContinuousMatch(this.minContinuousMatch);
    	lookup.loadPeptidesFromFile(peptideFile);
    	System.out.println("start searching");
    	long start = (new GregorianCalendar()).getTimeInMillis();
    	SimpleProbabilisticScorer scorer1 = (SimpleProbabilisticScorer)SpectrumUtil.getLPeakRankBaseScorer(training);
    	SimpleProbabilisticScorer filter = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedPeptideSingleScorer(mixtureTraining);
//    	SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(mixtureTraining);
    	SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(mixtureTraining);
    	((LinkedPeptidePeakScoreLearner)((MixtureSpectrumScorer)scorer2).getComp()).generalLinkedMode=true;
//    	List<Spectrum> specList = reader.readAllMS2Spectra();
//    	Iterator<Spectrum> iter = specList.iterator();
    	for(;iter.hasNext();){
    		Spectrum s = iter.next();
    		if(s.scanNumber < this.minScan || s.scanNumber > this.maxScan || s.charge >= 8){
    			continue;
    		}
    		//s.charge = (s.charge+1)*-1;
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
    		System.out.println(s.spectrumName + " has candidates: " + candidates.size());
    		List<Peptide> peps = LinkedPeakScoreLearner.generatePeptides(candidates);
    		List<Peptide> candPepWithMod = new ArrayList();
    		candPepWithMod.addAll(peps);
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 57, 1, 1));
    		//candPepWithMod.addAll(Peptide.insertPTM(peps, 28, 1, 1));
    		List<Peptide> linkedPeps = LinkedPeakScoreLearner.generateLinkedPeptides(candPepWithMod, s);
    		List<Spectrum> candidateSpectrum = LinkedPeakScoreLearner.generateSpectra3(linkedPeps, s);
    		//List<Spectrum> candidateSpectrum = LinkedPeakScoreLearner.generateSpectra2(candidates, s);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidateSpectrum, filter);
			if(targetPeptides.length == 2){
				int[] ranks = searcher.linkedRanks(s);
				System.out.println(s.spectrumName + " target peptides ranks " + ranks[0] + "\t" + ranks[1]);
			}
			Spectrum[] topSpectra = searcher.topSpectra(s, this.topFirstPassMatch);
			if(true){
			 	//continue;
			}
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
				double candParentMass = LookUpSpectrumLibX.getLinkedPartnerParentmass(p, s, CrossLinker.DSS);
				List<String> candidates2 = pLookup.getSpectrumByMass(candParentMass);
				//candidates2.addAll(pLookup.getSpectrumByMass(candParentMass-1.0));
				//candidates2.addAll(pLookup.getSpectrumByMass(candParentMass+1.0));
				for(int k = 0; k < candidates2.size(); k++){
					String pep2 = candidates2.get(k);
					Peptide p2 = new Peptide(candidates2.get(k) + ".2");
					double offset2 = LookUpSpectrumLibX.getLinkedOffSet(candidates2.get(k), s);
					//p2.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(o.getPeptide())-Mass.WATER});
					//o.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(candidates2.get(k))-Mass.WATER});
					int pos = pep2.indexOf(aaMatch);
					while(pos >= 0 ){//&& po != pep2.length()-1){
						Peptide copy = new Peptide(p2);
						copy.insertPTM(pos+1, offset2);
						//o.setCharge((short)2); //artifact from making linked peptide need to fix
						//System.out.println("peptide1: " + o + "\t"+ o.getCharge() + "\tPeptide2: "  + copy +"\t" + copy.getCharge());
						//TheoreticalSpectrum th = new TheoreticalSpectrum(o, copy, (short)s.charge, true);
						LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(o, copy, (short)s.charge);
						candidatePairs.add(th);
						pos = pep2.indexOf(aaMatch, pos+1);
					}
				}
				candidates2.add(p);
				Set<String> uniques = new HashSet<String>();
			}
			searcher = new SpectrumLibSearcher(candidatePairs, scorer2);
			searcher.setSingleScorer(scorer1);
			//searcher.topSpectrum(s);
			searcher.topLinkedSpectra(s, 10);
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
		mxdb.search(queryFile, peptideFile);
	}
	
	public static void main(String[] args){
		args[0] = "..\\mixture_linked\\MXDB_inputs2.txt";
		testMXDBSearch(args[0]);
	}
}
