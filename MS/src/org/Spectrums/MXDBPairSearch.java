package org.Spectrums;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Search for mixture and crosslinked peptides using pair peptide as input
 * @author bugouzhi
 *
 */
public class MXDBPairSearch {
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
		public MXDBPairSearch(){
			
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
//	    	SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(mixtureTraining);
	    	SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(mixtureTraining);
	    	//((LinkedPeptidePeakScoreLearner)((MixtureSpectrumScorer)scorer2).getComp()).generalLinkedMode=true;
//	    	List<Spectrum> specList = reader.readAllMS2Spectra();
//	    	Iterator<Spectrum> iter = specList.iterator();
	    	Mass.DSSLINKER_MASS = this.linkerMass;
	    	List<Peptide> pepList = getPeptidePair(peptideFile, new double[]{this.linkerMass, linkerMass+12.0753});	    	
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
	    		List<Spectrum> candidatePairs = new ArrayList<Spectrum>();
				lookup.setParentMassTolerance(3000.0);
				for(int i = 0; i < pepList.size(); i++){
					LinkedPeptide linkedPep = (LinkedPeptide)pepList.get(i);
					double parentmass = s.parentMass*s.charge - Mass.PROTON_MASS*(s.charge-1);
					//if(SpectrumUtil.checkMass(linkedPep.getParentmass(), parentmass, this.parentMassTolerance, SpectrumUtil.DALTON)){
					if(SpectrumUtil.checkMass(linkedPep.getParentmass(), parentmass, 12, SpectrumUtil.PPM)){
						LazyEvaluateLinkedSpectrum th = new LazyEvaluateLinkedSpectrum(linkedPep.peptides[0], linkedPep.peptides[1], (short)s.charge);
						candidatePairs.add(th);
					}
				}
				SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidatePairs, scorer2);
				searcher.matchTolerance = this.fragmentMassTolerance;
				searcher.setSingleScorer(scorer1);
				searcher.topSpectrum(s);
				searcher.topLinkedSpectra(s, 10, -1000);
				
	    	}
			System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		}
		
		public List<Peptide> getPeptidePair(String pepFile, double[] linkerMasses){
			List<String> pepPairs = Utils.FileIOUtils.createListFromFile(pepFile);
			List<Peptide> pepList = new ArrayList();
			for(int m = 0; m < linkerMasses.length; m++){
				CrossLinker linker = new CrossLinker(linkerMasses[m], 
						new int[]{CrossLinker.BUTCTERM}, new int[]{CrossLinker.BUTCTERM},
						new char[]{'K'}, new char[]{'K'});
				for(int i = 0; i < pepPairs.size(); i++){
					String[] tokens = pepPairs.get(i).split("\\s+");
					Peptide p1 = new Peptide(tokens[2], 1);
					Peptide p2 = new Peptide(tokens[3], 1);
					List<Peptide> peps = CandidateSpectrumLibFactory.getCrossLinkPeptides(p1, p2, linker, 1, 1);
					//List<LinkedPeptide> peps = linker.crossLinkPeptides(p1, p2, 1,1);
					pepList.addAll(peps);
				}
			}
			System.out.println("total pepetide pairs generated: " + pepList.size());
			return pepList;
		}
		
		public static void testMXDBSearch(String inFile){
			Map<String, String> arguments= Utils.FileIOUtils.createTableFromFile(inFile, 0, 1);
			MXDBPairSearch mxdb = new MXDBPairSearch();
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
			args[0] = "..\\mixture_linked\\MXDB_inputs21.txt";
			testMXDBSearch(args[0]);
		}
}

