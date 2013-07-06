package org.Spectrums;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * main interface to MixDB
 * @author Jian Wang
 *
 */
public class MixDB {
	public double windowWidth = 25;
	public int topPeaksKept = 10;
	//candidate filtering
	public int mode = 0;
	public double parentMassTolerance = 1;
	public double fragmentMassTolerance = 0.5;
	//input and training file
	public String queryFile;
	public String training = "..//mixture_linked//yeast_single_model_realannotated_win10_25.o";
	public String mixtureTraining = "";
	public int minScan = 0;
	public int maxScan = 0;
	
	
	public MixDB(){
		
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
    	SimpleProbabilisticScorer scorer1 = (SimpleProbabilisticScorer)SpectrumUtil.getRankBaseScorer(training);
    	SpectrumComparator scorer2 = SpectrumUtil.getMixtureScorer(mixtureTraining);
		List<Peak> pList = new ArrayList();
    	for(;iter.hasNext();){
    		Spectrum s = iter.next();
    		if(s.scanNumber < this.minScan || s.scanNumber > this.maxScan || s.charge >= 8){
    			continue;
    		}
    		//s.charge = (s.charge+1)*-1;
    		//s.windowFilterPeaks(this.topPeaksKept, this.windowWidth);
    		s.windowFilterPeaks(this.topPeaksKept, this.windowWidth);
    		s.removePrecursors(0.5);
    		s.computePeakRank();
    		System.out.println("Searching " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
			List<Peptide> chargedpeps = new ArrayList();
   			double lookupmass = s.parentMass*s.charge - Mass.PROTON_MASS*(s.charge-1) + Mass.PROTON_MASS;
   			List<PeptideLite> candidates = lookup.getCandidatePeptide(lookupmass, s.charge, pList);
   			List<Peptide> peps = LinkedPeakScoreLearner.generatePeptides(candidates);
   			for(int i = 0; i < peps.size(); i++){
   				Peptide p = peps.get(i);
   				for(int c = 2; c <= 3; c++){
   					if(Math.abs(s.parentMass - p.getParentmass()/c) < this.parentMassTolerance){
   						Peptide charged = new Peptide(p.getPeptide(), c);
   						chargedpeps.add(charged);
   					}
   				}
   			}
    		List<Spectrum> sList = new ArrayList();
    		for(Iterator<Peptide> it = peps.iterator(); it.hasNext();){
    			Peptide cand = it.next();
//    			LazyEvaluateLinkedSpectrum t = new LazyEvaluateLinkedSpectrum(((LinkedPeptide)cand).peptides[0], 
//						((LinkedPeptide)cand).peptides[1], cand.getCharge());
    			LazyEvaluatedSpectrum t = new LazyEvaluatedSpectrum(cand);
    			sList.add(t);
    		}
    		SpectrumLibSearcher searcher = new SpectrumLibSearcher(sList, scorer1, scorer2);
			searcher.setSingleScorer(scorer1);
			searcher.bestCandidates(s, 10);
    	}
 	}
	
	public static void testMixDB(String inFile){
		Map<String, String> arguments= Utils.FileIOUtils.createTableFromFile(inFile, 0, 1);
		MixDB mixdb = new MixDB();
		String[] prefixes = arguments.get("PrefixIons").split(",");
		String[] suffixes = arguments.get("SuffixIons").split(",");
		TheoreticalSpectrum.prefixIons = prefixes;
		TheoreticalSpectrum.suffixIons = suffixes;
		mixdb.mixtureTraining = arguments.get("MixtureTraining");
		mixdb.training = arguments.get("Training");
		mixdb.parentMassTolerance = Double.parseDouble(arguments.get("ParentMassTolerance"));
		mixdb.fragmentMassTolerance = Double.parseDouble(arguments.get("FragmentMassTolerance"));
		mixdb.topPeaksKept = Integer.parseInt(arguments.get("TopPeaksKept"));
		mixdb.windowWidth = Double.parseDouble(arguments.get("WindowWidth"));
		String queryFile = arguments.get("QueryFile");
		String peptideFile = arguments.get("PeptideFile");
		mixdb.mode = Integer.parseInt(arguments.get("ToleranceMode"));
		mixdb.minScan = Integer.parseInt(arguments.get("MinScan"));
		mixdb.maxScan = Integer.parseInt(arguments.get("MaxScan"));
		mixdb.search(queryFile, peptideFile);
	}
	
	public static void main(String[] args){
		args[0] = "..\\mixture_linked\\MXDBSumo_inputs.txt";
		testMixDB(args[0]);
	}
}
