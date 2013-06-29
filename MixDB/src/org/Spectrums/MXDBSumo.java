package org.Spectrums;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
	public double parentMassTolerance = 1;
	public double fragmentMassTolerance = 0.5;
	//input and training file
	public String queryFile;
	public String training = "..\\MSPLib\\Lib\\ecoli.msp";
	public String mixtureTraining = "";
	public int minScan = 0;
	public int maxScan = 0;
	
	
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
    	//SimpleProbabilisticScorer filter = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedPeptideSingleScorer(mixtureTraining);
    	SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(mixtureTraining);
    	//SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(mixtureTraining);
//   	 	((LinkedPeptidePeakScoreLearner)((MixtureSpectrumScorer)scorer2).getComp()).generalLinkedMode=true;
		List<Peak> pList = new ArrayList();
		//creat sumo
		List<Peptide> peptides2 = new ArrayList();
		peptides2.add(new Peptide("Z", 1));
		peptides2.add(new Peptide("Q-17.027QQTGG", 1));
//		peptides2.add(new Peptide("Q-17.027EQTGG", 1));
//		peptides2.add(new Peptide("Q-17.027QPTGG", 1));
		peptides2.add(new Peptide("QQQTGG", 1));
//		peptides2.add(new Peptide("QEQTGG", 1));
//		peptides2.add(new Peptide("QQPTGG", 1));
//		peptides2.add(new Peptide("EMEDEDTIDVFQQPTGG", 1));
//		peptides2.add(new Peptide("EMEDEDTIDVFQQQTGG", 1));
//		peptides2.add(new Peptide("GMEEEDVIEVYQEQTGG", 1));
//		peptides2.add(new Peptide("EDEDTIDVFQQPTGG", 1));
//		peptides2.add(new Peptide("EDEDTIDVFQQQTGG", 1));
//		peptides2.add(new Peptide("EEEDVIEVYQEQTGG", 1));		
		Mass.DSSLINKER_MASS = -1*Mass.WATER;
		CrossLinker SUMO = new CrossLinker(Mass.WATER*-1, new int[]{CrossLinker.ANYPOSITION}, new int[]{CrossLinker.CTERM, CrossLinker.ANYPOSITION}, new char[]{'K'}, new char[]{'G','Z'});
    	for(;iter.hasNext();){
    		Spectrum s = iter.next();
    		if(s.scanNumber < this.minScan || s.scanNumber > this.maxScan || s.charge >= 8){
    			continue;
    		}
    		//s.charge = (s.charge+1)*-1;
    		//s.windowFilterPeaks(this.topPeaksKept, this.windowWidth);
    		s.windowFilterPeaks(8, 25);
    		s.removePrecursors(0.5);
    		s.computePeakRank();
    		System.out.println("Searching " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
			List<Peptide> linkedpeps = new ArrayList();
    		for(int i = 0; i < peptides2.size(); i++){
    			Peptide p2 = peptides2.get(i);
    			double lookupmass = s.parentMass*s.charge - Mass.PROTON_MASS*(s.charge-1) - p2.getParentmass() - SUMO.getLinkerMassOffSet() + Mass.PROTON_MASS;
    			List<PeptideLite> candidates = lookup.getCandidatePeptide(lookupmass, s.charge, pList);
    			List<Peptide> peps = LinkedPeakScoreLearner.generatePeptides(candidates);
    			//double lookupmass2 = s.parentMass*s.charge - Mass.PROTON_MASS*(s.charge-1) - p2.getParentmass() - SUMO.getLinkerMassOffSet() + Mass.PROTON_MASS - 42.010565;
    			//List<PeptideLite> candidates2 = lookup.getCandidatePeptide(lookupmass2, s.charge, pList);
    			//List<Peptide> peps2 = LinkedPeakScoreLearner.generatePeptides(candidates2);
    			//peps2 = Peptide.insertPTM(peps2, 43.025, 1, 1);
    			//peps2 = Peptide.insertPTM(peps2, 42.010565	, 1, 1);
    			//peps.addAll(peps2);
    			linkedpeps.addAll(CandidateSpectrumLibFactory.getCrossLinkPeptides(peps, p2, SUMO, s.charge, s.charge));
    		}
    		List<Spectrum> sList = new ArrayList();
    		for(Iterator it = linkedpeps.iterator(); it.hasNext();){
    			LinkedPeptide cand = (LinkedPeptide)it.next();
    			LazyEvaluateLinkedSpectrum t = new LazyEvaluateLinkedSpectrum(((LinkedPeptide)cand).peptides[0], 
						((LinkedPeptide)cand).peptides[1], cand.getCharge());
    			sList.add(t);
    		}
    		SpectrumLibSearcher searcher = new SpectrumLibSearcher(sList, scorer2);
			searcher.setSingleScorer(scorer1);
			searcher.topLinkedSpectra(s, 10);
    	}
 	}
	
	public static void testMXDBSUMO(String inFile){
			Map<String, String> arguments= Utils.FileIOUtils.createTableFromFile(inFile, 0, 1);
		MXDBSumo mxdb = new MXDBSumo();
		String[] prefixes = arguments.get("PrefixIons").split(",");
		String[] suffixes = arguments.get("SuffixIons").split(",");
		TheoreticalSpectrum.prefixIons = prefixes;
		TheoreticalSpectrum.suffixIons = suffixes;
		mxdb.mixtureTraining = arguments.get("MixtureTraining");
		mxdb.training = arguments.get("Training");
		mxdb.parentMassTolerance = Double.parseDouble(arguments.get("ParentMassTolerance"));
		mxdb.fragmentMassTolerance = Double.parseDouble(arguments.get("FragmentMassTolerance"));
		mxdb.topPeaksKept = Integer.parseInt(arguments.get("TopPeaksKept"));
		mxdb.windowWidth = Double.parseDouble(arguments.get("WindowWidth"));
		String queryFile = arguments.get("QueryFile");
		String peptideFile = arguments.get("PeptideFile");
		mxdb.mode = Integer.parseInt(arguments.get("ToleranceMode"));
		mxdb.minScan = Integer.parseInt(arguments.get("MinScan"));
		mxdb.maxScan = Integer.parseInt(arguments.get("MaxScan"));
		mxdb.search(queryFile, peptideFile);
	}
	
	public static void main(String[] args){
		args[0] = "..\\mixture_linked\\MXDBSumo_inputs.txt";
		testMXDBSUMO(args[0]);
	}
}
