package org.Spectrums;
import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar; 
import java.util.Iterator;
import java.util.List;

public class testSpectrumLibSearch {
	public static void testRank(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		String file = "..\\mixture_linked\\Ecoli_allpeptides.txt";
		lib1.removeModSpectra();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		List<Spectrum>  specList = lib1.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < 50; i++){
			Spectrum s = specList.get(i);
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 5, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), scorer);
			System.out.println("Spectrum: " + s.peptide + " rank " + searcher.rank(s) + " in spectrumLib of size: " + lib.getSpectrumList().size());			
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void testLinkedRank(){
		String spectrumFile = ".\\mixture_linked\\spectrums_raw.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
		String file = "..\\mixture_linked\\longrunInpsect_plus_allEcoli_peptides.txt";
		lib1.removeModSpectra();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		List<Spectrum>  specList = lib1.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < 30; i++){
			Spectrum s = specList.get(i);
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 5, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), scorer);
			System.out.println("Spectrum: " + s.peptide + " rank " + searcher.rank(s) + " in spectrumLib of size: " + lib.getSpectrumList().size());			
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void removeSUMOPeaks(Spectrum s){
		Peptide peptide = new Peptide("QQQTGG", 1);
		List<Peptide> pepList = new ArrayList<Peptide>();
		pepList.add(peptide);
		List<Peptide> linkedPeps = LinkedPeakScoreLearner.generateLinkedPeptides(pepList, s);
		Peptide linkPep = linkedPeps.get(1);
		TheoreticalSpectrum t = new TheoreticalSpectrum(linkPep, s.charge);
		SpectrumUtil.removeAnnotatedPeaks(s, t, 0.5);
	}
	
	public static void testTopMatch(String spectrumFile, String spectrumFile2, String file, double tolerance){
		//spectrumFile = "..\\MSPLib\\Lib\\yeast.msp";
		spectrumFile2 = "../mixture_linked/linked_peptide_library/sumo_lib/human_sumo/20110527_Soufib_Soufib_Pscan_SUMO_H2O2_CID.mzXML";
		tolerance = 3.0;
		file = "..\\mixture_linked\\testpeptides.txt";
		//SpectrumLib lib2 = new SpectrumLib("..\\MSPLib\\Lib\\Ecoli.msp", "MSP");
		//lib2.removeModSpectra();
		//String probFile = ".\\data\\IonsScore.txt";
		//String noiseModel = ".\\data\\NoiseModel.txt";
		//SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		//SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		//List<Spectrum>  specList = SpectrumUtil.getRandomSpectrum(lib2, 100);
		//Iterator<Spectrum> iterator = lib2.getAllSpectrums().iterator();
		MZXMLReader iterator = new MZXMLReader(spectrumFile2);
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		factory.setMinCharge(2);
		factory.setMaxCharge(4);
		factory.matchCharge = true;
		factory.indexPeptideByParentMass(0.5);
		factory.insertPTM(43.025, 1, 1);
		//factory.insertPTM(28, 1, 1);
		factory.insertPTM(551.230, new char[]{'K'}, 2);
		//factory.insertPTM(582.239, new char[]{'K'}, 1);
		//factory.insertPTM(-57, new char[]{'C'}, 1);
		//factory.insertPTM(42.010565, 1, 2);
		LabelledPeakFactory.clearFactory();
		PeakComparator pComp = RankBaseScoreLearner.loadComparator("../mixture_linked/yeast_single_peptide_model.o");
		SpectrumComparator scorer1 = new SimpleProbabilisticScorer(pComp);
		//lib2 = null;
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(;iterator.hasNext();){
			Spectrum s = iterator.next();
			//s.windowFilterPeaks(30, 50);
			if(s.scanNumber != 3232){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.removePrecursors(0.5);
			//removeSUMOPeaks(s);
			System.out.println(s.spectrumName + "\t" + s.charge + "\tobserved spectrum has peaks: " + s.getPeak().size());
			s.computePeakRank();
			if(s.getPeaks().size() < 10){
				continue;
			}
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, tolerance, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), scorer1);
			searcher.bestSpectra(s, 10);
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMixtureRank(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		String file = "..\\mixture_linked\\Yeast_allpeptides_plusSpecLib.txt";
		lib1.removeModSpectra();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		//String mixturefile = "..\\mixture_linked\\yeast_mixture.name";
		String mixturefile = "..\\mixture_linked\\exp_sim_mixtures.id";
		String spectrumFile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
//		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
//		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		SpectrumLib mixlib = lib2.createMix(mixturefile, 1000, 1, 1.0, 0.000000001, 1.0, 2000, false);
		//SpectrumLib mixlib = lib1.createRandomMix(1000, 1, 1, 0.0005, 1.0, 2, false);
		//List<Spectrum>  specList = SpectrumUtil.getRandomSpectrum(lib1, 1000);
		lib1 = null;
		List<Spectrum>  specList = mixlib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			//s.peptide = s.peptide + " & " + s.peptide;
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), filter, filter);
			int[] ranks = searcher.ranks(s);
			System.out.println("Spectrum: " + s.peptide + " ranks " +  ranks[0] + " & " + ranks[1]+ " in spectrumLib of size: " + lib.getSpectrumList().size());
			//searcher.bestPair(s);
			//searcher.bestCandidates(s, 10);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMixtureRankBaseScorer(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		String mixFile = "..\\mixture_linked\\linked_peptide_spectra2.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		//String file = "..\\mixture_linked\\Ecoli_allpeptides.txt";
		String file = "..\\mixture_linked\\Ecoli_allpeptides_plusLinkedpeptides.txt";
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		//String mixturefile = "..\\mixture_linked\\yeast_mixture.name"; 
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		
//		SpectrumIonRankLearner learner = new SpectrumIonRankLearner(lib1);
//		PeakComparator peakscorer2 = learner.createComparatorSet();
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer2);
		//RankBaseMixtureSpectrumScorer scorer = new RankBaseMixtureSpectrumScorer(peakscorer2);
		
		//SpectrumLib mixlib = lib1.createMix(mixturefile, 1000, 1, 0.0, 0.00001, 1.0, 2000, true);
		//SpectrumLib mixlib = lib1.createRandomMix(1000, 1, 0.3, 0.0005, 1.0, 2, false);
		//SpectrumLib mixlib = lib1.createRandomMix(1000, 0.1, 0.0005, 1.0, 2, false);
		//List<Spectrum>  specList = SpectrumUtil.getRandomSpectrum(lib1, 1000);
		lib1 = null;
		//List<Spectrum>  specList = lib1.getSpectrumList();
		SpectrumLib mixlib = new SpectrumLib(mixFile, "MGF");
		//mixlib.windowFilterPeaks(7, 25);
		//mixlib.computeRank();
		List<Spectrum>  specList = mixlib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			//s.peptide = s.peptide + " & " + s.peptide;
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 1.5, false);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), filter, scorer);
			//int[] ranks = searcher.ranks(s);
			//System.out.println("Spectrum: " + s.peptide + " ranks " +  ranks[0] + " & " + ranks[1]+ " in spectrumLib of size: " + lib.getSpectrumList().size());
			searcher.bestCandidates(s, 10);
			//searcher.bestSpectrum(s);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMixtureSearchMethod(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		//String file = "..\\mixture_linked\\longrunInpsect_plus_allEcoli_peptides.txt";
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusDecoy.txt";
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		//SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer2);
		RankBaseMixtureSpectrumScorer scorer = new RankBaseMixtureSpectrumScorer(peakscorer2);
		
		//String mixlibFile = "..\\mixture_compressed\\new80min.mgf";
		String mixlibFile = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
		//SpectrumLib mixlib = new SpectrumLib(mixlibFile, "MGF");
		Iterator<Spectrum> specIterator = new LargeSpectrumLibIterator<Spectrum>(mixlibFile);
		lib1 = null;
		//mixlib.windowFilterPeaks(7, 25);
		//mixlib.computeRank();
		//List<Spectrum>  specList = mixlib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		
		for(;specIterator.hasNext();){
			Spectrum s = specIterator.next();
//			String[] tokens = s.spectrumName.split("[_.]");
//			int num = Integer.parseInt(tokens[1]);
//			if(num < 80 || num > 400){
//				//continue;
//			}
			s.windowFilterPeaks(10, 25);
			if(s.getPeak().size() < 30){
				continue;
			}
			s.computePeakRank();	
			//s.peptide = s.peptide + " & " + s.peptide;
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			System.out.println("After mFilter we have " + lib.getAllSpectrums().size());
//			LookUpSpectrumLib filteredLib = new LookUpSpectrumLib(lib.getSpectrumList());
//			List<Peak> queryPeaks = s.getTopPeaks(25);
//			//List<Peak> queryPeaks = s.topWindowPeaks(6, 200);
//			System.out.println("Query peaks has: " + queryPeaks.size());
//			List<Spectrum> candidates = filteredLib.getSpectrumByPeaks(queryPeaks, 4);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getAllSpectrums(), filter, scorer);
			//SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getSpectrumList(), scorer);

			//int[] ranks = searcher.ranks(s);
			//System.out.println("Spectrum: " + s.peptide + " ranks " +  ranks[0] + " & " + ranks[1]+ " in spectrumLib of size: " + lib.getSpectrumList().size());
			searcher.bestCandidates(s, 10);
			//searcher.topSpectrum(s);
			//searcher.bestSpectrum(s);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testLookUpTableSearch(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		String spectrumFile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		//String file = "..\\mixture_linked\\Yeast_allPeptides_plusDecoy.txt";
		String file = "..\\mixture_linked\\Yeast_allPeptides_plusSpecLib_plusDecoy.txt";
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		String mixturefile = "..\\mixture_linked\\ecoli_mixture2.name";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		
		//SpectrumIonRankLearner learner = new SpectrumIonRankLearner(lib1);
		//PeakComparator peakscorer2 = learner.createComparatorSet();
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib1);
		peakscorer2.getIonsCount();
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer2);
		//RankBaseMixtureSpectrumScorer scorer = new RankBaseMixtureSpectrumScorer(peakscorer2);
		
		//SpectrumLib mixlib = lib1.createMix(mixturefile, 1000, 1, 0.0, 0.00001, 1.0, 2000, true);
		//SpectrumLib mixlib = lib1.createRandomMix(2000, 1, 0.0, 0.0005, 1.0, 2, false);
		//SpectrumLib mixlib = lib1.createRandomMix(1000, 0.1, 0.0005, 1.0, 2, true);
		//List<Spectrum>  specList = SpectrumUtil.getRandomSpectrum(lib1, 1000);
		//lib1 = null;
		//List<Spectrum>  specList = lib1.getSpectrumList();
		//SpectrumLib mixlib = new SpectrumLib(mixFile, "MGF");
		//mixlib.windowFilterPeaks(7, 25);
		//mixlib.computeRank();
		//List<Spectrum>  specList = mixlib.getSpectrumList();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		long start = (new GregorianCalendar()).getTimeInMillis();
		Iterator<Spectrum> specIterator = new LargeSpectrumLibIterator<Spectrum>(spectrumFile2);
		int i =0;
		for(;specIterator.hasNext();){
			Spectrum s = specIterator.next();
			if(s.getPeak().size() < 30){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			SpectrumLib lib = factory.createCandidateSpectrumLibX(s, 3, false);
			System.out.println("After filter one we have: " + lib.getAllSpectrums().size() + " candidates");
			LookUpSpectrumLib filteredLib = new LookUpSpectrumLib(lib.getSpectrumList());
			List<Peak> queryPeaks = s.getTopPeaks(25);
			//List<Peak> queryPeaks = s.topWindowPeaks(6, 200);
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<Spectrum> candidates = filteredLib.getSpectrumByPeaks(queryPeaks, 4);
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidates, scorer, scorer);
			//int[] ranks = searcher.ranks(s);			
			//System.out.println("Spectrum: " + s.peptide + " ranks " +  ranks[0] + " & " + ranks[1]+ " in spectrumLib of size: " + lib.getSpectrumList().size());
			searcher.topSpectrum(s);
			i++;
			if(i > 200){
				break;
			}
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void testSumoSearch(String queryFile, String peptidePattern, String peptideFile, String training, String training2, int position1, int position2, int minScan, int maxScan){
		CombinatoryPeptides combPeps = new CombinatoryPeptides(peptidePattern); 
		System.out.println("pattern is: peptidePattern");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("TALH[ARDEHLKMFPSTYV]K[ARDEHLKMFPSTYV]S[ARDEHLKMFPSTYV]TFR");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("A[FILVYW]K[ARDEHLKMFPSTYV][DE]T[ARDEHLKMFPSTYV]FRAK"); 
		List<String> peptides = combPeps.generateAllPeptides();
		peptides.clear();
		List<String> ecoPep = Utils.FileIOUtils.createListFromFile(peptideFile);
		peptides.addAll(ecoPep);
		CombinatoryPeptides combPeps2 = new CombinatoryPeptides("Z");
		List<String> peptides2 = combPeps2.generateAllPeptides();
		peptides2.add("Q-17.027QQTGG");
		peptides2.add("QQQTGG");
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);	
		factory.indexPeptideByParentMass(1.0);
		//factory.insertPTM(42, 1, 1);
		CandidateSpectrumLibFactory factory2 = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides2);
		factory2.setMinCharge(1);
		factory2.setMaxCharge(1);	
		factory2.indexPeptideByParentMass(1.0);
		CrossLinker SUMO = new CrossLinker(Mass.WATER*-1, new int[]{CrossLinker.ANYPOSITION}, new int[]{CrossLinker.CTERM, CrossLinker.ANYPOSITION}, new char[]{'K'}, new char[]{'G','Z'});
		System.out.println("start crosslinking");
		factory.crossLinkAllPeptides(factory2, 2, 4, SUMO);
		System.out.println("done crosslinking");
		//factory.crossLinkAllPeptides(2,3);
		Iterator<Spectrum> it = null;
		if(queryFile.contains("mzXML")){
				MZXMLReader reader = new MZXMLReader(queryFile);
				it = reader;
		}
		if(queryFile.contains("mgf")){
				it = new LargeSpectrumLibIterator(queryFile);
		}
		//List<Spectrum> specList = reader.readAllMS2Spectra();
		//SpectrumLib queries = new SpectrumLib(spectrumLibFile, "MGF");
		//List<Spectrum> specList = queries.getAllSpectrums();
		long start = (new GregorianCalendar()).getTimeInMillis();
		SpectrumComparator scorer1 = SpectrumUtil.getLPeakRankBaseScorer(training);
		SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(training2);
		factory.setMatchCharge(true);
		for(; it.hasNext();){
			Spectrum s = it.next();
			//s.windowFilterPeaks(30, 50);
			//System.out.println("scan number is: " + s.scanNumber);
			if(s.parentMass < 400 || s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.removePrecursors(0.5);
			s.computePeakRank();
			List<Peptide> candidates = factory.getCandidateByMass(s.parentMass, s.charge, 3);
			//List<Peptide> candidates2 = factory.getCandidateByMass(s.parentMass-(Mass.C13-Mass.C12)/s.charge, s.charge, 0.05);
			//candidates.addAll(candidates2);
			System.out.println(s.spectrumName + "\t" + s.parentMass + "\t" +  s.charge + "\t" + " observed spectrum has peaks: " + s.getPeak().size() + " number of candidate peptides: " + candidates.size());
			List<Spectrum> slist = new ArrayList<Spectrum>();
			for(int j = 0; j < candidates.size(); j++){
				slist.add(CombinatoryPeptides.generateLinkedTheoreticalSpectrum((LinkedPeptide)candidates.get(j)));
			}
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(slist, scorer2);
			searcher.setSingleScorer(scorer1);
			searcher.topLinkedSpectra(s, 10);
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");		
	}
	
	public static CandidateSpectrumLibFactory getSUMOFactory(String peptideFile){
		List<String> peptides = Utils.FileIOUtils.createListFromFile(peptideFile);
		CombinatoryPeptides combPeps2 = new CombinatoryPeptides("Z");
		List<String> peptides2 = combPeps2.generateAllPeptides();
		peptides2.add("Q-17.027QQTGG");
		peptides2.add("QQQTGG");
		peptides2.add("QQTGG");
		peptides2.add("QTGG");
		peptides2.add("TGG");
		peptides2.add("GG");
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setPeptideFile(peptideFile);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);	
		factory.indexPeptideByParentMass(1.0);
		factory.insertPTM(42.010565, 1, 1);
		//factory.insertPTM(0.984, new char[]{'N', 'Q'}, 1);
		//factory.insertPTM(6.01, new char[]{'V', 'P'}, 2);
		//factory.insertPTM(7.01, new char[]{'L'}, 2);
		factory.insertPTM(-57, new char[]{'C'}, 2);
		CandidateSpectrumLibFactory factory2 = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides2);
		factory2.setMinCharge(1);
		factory2.setMaxCharge(1);	
		factory2.indexPeptideByParentMass(1.0);
		CrossLinker SUMO = new CrossLinker(Mass.WATER*-1, new int[]{CrossLinker.ANYPOSITION}, new int[]{CrossLinker.CTERM, CrossLinker.ANYPOSITION}, new char[]{'K'}, new char[]{'G','Z'});
		System.out.println("start crosslinking");
		factory.crossLinkAllPeptides(factory2, 2, 4, SUMO);
		System.out.println("done crosslinking");
		return factory;
	}
	
	public static CandidateSpectrumLibFactory getCrossLinkFactory(String peptideFile){
		List<String> peptides = Utils.FileIOUtils.createListFromFile(peptideFile);
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setPeptideFile(peptideFile);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);	
		factory.indexPeptideByParentMass(1.0);
		factory.insertPTM(42.010565, 1, 1);
		//factory.insertPTM(28, 1, 1);
		//factory.insertPTM(156, new char[]{'K'}, 1);
		//factory.insertPTM(1, new char[]{'Q'}, 1);
		Mass.DSSLINKER_MASS = -116.0430;//136.0995;//138.0680;
		CrossLinker CROSSLINKER = new CrossLinker(Mass.DSSLINKER_MASS, new int[]{CrossLinker.ANYPOSITION}, new int[]{CrossLinker.ANYPOSITION}, new char[]{'C'}, new char[]{'C'});
		System.out.println("start crosslinking");
		factory.crossLinkAllPeptides(factory, 3, 5, CROSSLINKER);
		//factory.internalLinkPeptides(2, 7, CROSSLINKER);
		System.out.println("done crosslinking");
		return factory;
	}
	
	public static void runSUMOSearch(String queryFile, String peptideFile, String scoreFile1, String scoreFile2, double precursor, int minScan, int maxScan){
		TheoreticalCandidateSpectrumFactory factory = null;
		if(peptideFile.contains(".map")){
			factory = new TheoreticalCandidateSpectrumFactory(peptideFile);
		}else{
			factory = new TheoreticalCandidateSpectrumFactory(getSUMOFactory(peptideFile));
		}
		Iterator<Spectrum> it = null;
		if(queryFile.contains("mzXML")){
				SortedMZXMLReader reader = new SortedMZXMLReader(queryFile);
				it = reader;
		}
		if(queryFile.contains("mgf")){
				it = new LargeSpectrumLibIterator(queryFile);
		}
		//List<Spectrum> specList = reader.readAllMS2Spectra();
		//SpectrumLib queries = new SpectrumLib(spectrumLibFile, "MGF");
		//List<Spectrum> specList = queries.getAllSpectrums();
		long start = (new GregorianCalendar()).getTimeInMillis();
		SpectrumComparator scorer1 = SpectrumUtil.getLPeakRankBaseScorer(scoreFile1);
		//SimpleProbabilisticScorer scorer1 = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedPeptideSingleScorer(scoreFile2);
		SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(scoreFile2);
		//SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(scoreFile2);
		factory.setMatchCharge(true);
		for(; it.hasNext();){
			Spectrum s = it.next();
			//s.windowFilterPeaks(30, 50);
			//System.out.println("scan number is: " + s.scanNumber);
			if(s.parentMass < 400 || s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.removePrecursors(0.5);
			s.computePeakRank();
			List candidates = new ArrayList();
			candidates.addAll(factory.getSUMOCandidates(s, precursor));
			//List<Peptide> candidates2 = factory.getCandidateByMass(s.parentMass-(Mass.C13-Mass.C12)/s.charge, s.charge, 0.05);
			System.out.println(s.spectrumName + "\t" + s.parentMass + "\t" +  s.charge + "\t" + " observed spectrum has peaks: " + s.getPeak().size() + " number of candidate peptides: " + candidates.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidates, scorer2);
			searcher.setSingleScorer(scorer1);
			searcher.topLinkedSpectra(s, 10);
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void runBruteForaceCrossLinkSearch(String queryFile, String peptideFile, String scoreFile1, String scoreFile2, int minScan, int maxScan){
		TheoreticalCandidateSpectrumFactory factory = null;
		if(peptideFile.contains(".map")){
			factory = new TheoreticalCandidateSpectrumFactory(peptideFile);
		}else{
			factory = new TheoreticalCandidateSpectrumFactory(getCrossLinkFactory(peptideFile));
		}
		Iterator<Spectrum> it = null;
		if(queryFile.toLowerCase().contains("mzxml")){
				SortedMZXMLReader reader = new SortedMZXMLReader(queryFile);
				it = reader;
		}
		if(queryFile.contains("mgf")){
				it = new LargeSpectrumLibIterator(queryFile);
		}
		//List<Spectrum> specList = reader.readAllMS2Spectra();
		//SpectrumLib queries = new SpectrumLib(spectrumLibFile, "MGF");
		//List<Spectrum> specList = queries.getAllSpectrums();
		long start = (new GregorianCalendar()).getTimeInMillis();
		SpectrumComparator scorer1 = SpectrumUtil.getLPeakRankBaseScorer(scoreFile1);
		//SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer(scoreFile2);
		//((LinkedPeptidePeakScoreLearner)((MixtureSpectrumScorer)scorer2).getComp()).generalLinkedMode=true;
		SpectrumComparator scorer2 = SpectrumUtil.getLMixtureScorer(scoreFile2);
		factory.setMatchCharge(true);
		for(; it.hasNext();){
			Spectrum s = it.next();
			//s.windowFilterPeaks(30, 50);
			//System.out.println("scan number is: " + s.scanNumber);
			if(s.charge > 8 || s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}
			s.windowFilterPeaks(8, 25);
			s.removePrecursors(0.5);
			s.computePeakRank();
			List candidates = new ArrayList();
			candidates.addAll(factory.getSUMOCandidates(s, 1));
			//List<Peptide> candidates2 = factory.getCandidateByMass(s.parentMass-(Mass.C13-Mass.C12)/s.charge, s.charge, 0.05);
			System.out.println(s.spectrumName + "\t" + s.parentMass + "\t" +  s.charge + "\t" + " observed spectrum has peaks: " + s.getPeak().size() + " number of candidate peptides: " + candidates.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidates, scorer2);
			searcher.setSingleScorer(scorer1);
			searcher.topLinkedSpectra(s, 10);
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void main(String[] args){
		//testRank();
		//testMixtureRank();
		//testMixtureRankBaseScorer();	
		//testMixtureSearchMethod();	
		//testTopMatch(args[0], args[1], args[2], Double.parseDouble(args[3]));
		//testLookUpTableSearch();
		//args[0] = "..\\mixture_linked\\lib_sumo1_spectra_long_single_peptide_search_with_pyroQ.mgf";
		//args[0] = "../mixture_linked/linked_peptide_library/toni/110520_Disulfide/tk110517_Nuno_DMSO_peptides3.mzXML";
		args[0] = "..\\mixture_linked/linked_peptide_library/sumo_lib/mcl1/MCL1_Bo_Light_200fm_40%Collision.mzXML";
//		//args[1] = "KA[ARDEHLKMFPSTYV]D[ARDEHLKMFPSTYV]ES[ARDEHLKMFPSTYV]LRAK";
//		//args[1] = "A[FILVYW]K[ARDEHLKMFPSTYV][DE]T[ARDEHLKMFPSTYV]FRAK";
		args[1] = "..\\mixture_linked\\testpeptides.txt";
		//args[1] = "..\\mixture_linked\\database/lib_sumo_peptides_plusEcoli_allPeptides.txt_processed.map";
		args[2] = "..\\mixture_linked\\yeast_linked_single_peptide_model.o";
		args[3] = "..\\mixture_linked\\lib_sumo_linked_score_model1.o";
		//args[3] = "..\\mixture_linked\\linked_mixture_model_alpha0.5.o";
		args[4] = "0.05";
		args[5] = "788";
		args[6] = "788";
//		//testSumoSearch(args[0], args[1], args[2], args[3], args[4], Integer.parseInt(args[5]), Integer.parseInt(args[6]), Integer.parseInt(args[7]), Integer.parseInt(args[8]));
		runSUMOSearch(args[0], args[1], args[2], args[3], Double.parseDouble(args[4]), Integer.parseInt(args[5]), Integer.parseInt(args[6]));
//		runBruteForaceCrossLinkSearch(args[0], args[1], args[2], args[3], Integer.parseInt(args[5]), Integer.parseInt(args[6]));
	}
}

