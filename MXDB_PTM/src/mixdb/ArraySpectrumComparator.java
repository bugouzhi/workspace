package mixdb;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;

import org.Spectrums.ArrayUtils;
import org.Spectrums.LinkedPeptide;
import org.Spectrums.LinkedPeptidePeakScoreLearner;
import org.Spectrums.MZXMLReader;
import org.Spectrums.MixturePeakScoreLearner;
import org.Spectrums.PeakComparator;
import org.Spectrums.Peptide;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.SimpleProbabilisticScorer;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumComparator;
import org.Spectrums.SpectrumLib;
import org.Spectrums.TheoreticalSpectrum;

public class ArraySpectrumComparator implements SpectrumComparator{
	private double massTolerance = 0.5;
	private double[][] scoreTable;    //note to self: should we consider having score table directly here or an adaptor object instead??
	private double[][] errorTable;    //maybe it is still better to use score table directly, since adaptor don't need if we completely
	private PeakComparator pComp;     //migrate to the new framework
	private boolean DEBUG = false;
	private IonTypeMapper ionMap;
	private double[] massErrorInterval = {-0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55};
	
	
	public ArraySpectrumComparator(double tolerance, PeakComparator pComp){
		this.massTolerance = tolerance;
		this.pComp = pComp;
	}
	@Override
	public double compare(Spectrum s1, Spectrum s2) {
		if(!(s1 instanceof ArrayTheoreticalSpectrum)){
			throw new IllegalArgumentException("spectrum 1 being compared must be of type ArrayTheoreticalSpectrum");
		}
		
		if(!(s2 instanceof ArraySpectrum)){
			throw new IllegalArgumentException("spectrum 2 being compared must be of type ArraySpectrum");
		}
		double score = 0.0, matchesScore = 0.0, unMatchesScore = 0.0;
		ArraySpectrum array1 = (ArraySpectrum)s1;
		ArraySpectrum array2 = (ArraySpectrum)s2;
		double[][] massIntList1 = array1.getMassIntensityList();
		double[][] massIntList2 = array2.getMassIntensityList();
		boolean[] isMatched = new boolean[massIntList1[ArraySpectrum.MASS].length];
		int matchedCount = 0, matchedCount2=0, unMatchedCount = 0;
		int i = 0, j = 0;
		while(j < massIntList2[0].length && i < massIntList1[0].length){
			//System.out.println(i + "\t" + j);
			double mass2 = massIntList2[ArraySpectrum.MASS][j];
			double mass1 = massIntList1[ArraySpectrum.MASS][i];
			if(mass2 - mass1 > this.massTolerance){
				i++;
				continue;
			}
			if(mass1 - mass2 > this.massTolerance){
				j++;
				continue;
			}
			double peakMaxScore = -10000.0;
			double mass = mass1;
//			for(int p = i; mass - mass2 < this.massTolerance && p < massIntList1[0].length; 
//				p++, mass = massIntList1[ArraySpectrum.MASS][p]){
			int p = i;
			while(p < massIntList1[0].length){
				mass = massIntList1[ArraySpectrum.MASS][p];
				if(mass - mass2 > this.massTolerance){
					break;
				}
				isMatched[p]=true;
				double currScore = 0.0;
				double peakscore = 	getPeakScore((int)massIntList1[ArraySpectrum.INTENSITY][p],
						(int)massIntList2[ArraySpectrum.INTENSITY][j]);
				double errorScore = getErrorScore(massIntList1[ArraySpectrum.MASS][p]-massIntList2[ArraySpectrum.MASS][j],
						(int)massIntList2[ArraySpectrum.INTENSITY][j]);
				currScore = peakscore + errorScore;
				peakMaxScore = peakMaxScore > currScore ? peakMaxScore : currScore;
				if(DEBUG){
					System.out.println("peak: " + massIntList2[0][j] + "\t" + massIntList2[1][j] + "\t~\t" 
							+  this.ionMap.getIonType((int)massIntList1[1][p])+ "\tcurrent score: " 
							+ peakscore + "\t"  + errorScore + "\t" + currScore);
				}
				p++;
			}
			
			if(peakMaxScore != -10000){
				matchedCount++;
			}
			matchesScore += peakMaxScore;
			if(DEBUG){
				System.out.println("current max: " + peakMaxScore);
			}
			j++;
			continue;
		}
		//take care of unmatched theoretical fragments
		for(int k = 0; k < isMatched.length; k++){
			if(!isMatched[k]){
				unMatchesScore += getPeakScore((int)massIntList1[ArraySpectrum.INTENSITY][k], 0);  //zero means unmatched peaks
				unMatchedCount++;
			}else{
				matchedCount2++;
			}
		}
		score = matchesScore + unMatchesScore;
		if(DEBUG){
			System.out.println("score: " + score + "\t" + matchesScore + "\t" + unMatchesScore);
			System.out.println("matched count: " + matchedCount + "(" + matchedCount2 + ")" + "\t" + unMatchedCount);
		}
		return score;
	}
	
	public double getPeakScore(int type, int rank){
		//System.out.println("type: " + type + "\t" + "ranks: " + rank);
		return this.scoreTable[type][rank];
	}
	
	public double getErrorScore(double error, int rank){
		if(DEBUG){
			System.out.println("rank is: " + rank + " error is: " + error);
		}
		int errorInd = ArrayUtils.getIntervalIndex(error, this.massErrorInterval);
		return this.errorTable[rank][errorInd];
	}
	
	public double[][] getScoreTable() {
		return scoreTable;
	}
	public void setScoreTable(double[][] scoreTable) {
		this.scoreTable = scoreTable;
	}
	public double[][] getErrorTable() {
		return errorTable;
	}
	public void setErrorTable(double[][] errorTable) {
		this.errorTable = errorTable;
	}
	
	public static ArraySpectrumComparator loadStandardComparator(String path){
		RankBaseScoreLearner peakscorer = RankBaseScoreLearner.loadComparatorLocal(path);
		RankBaseScoreAdapter adaptor = new RankBaseScoreAdapter(TheoreticalSpectrumFactory.standardTypeMap, 
				TheoreticalSpectrumFactory.standardIonMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		arryScorer.scoreTable = adaptor.getTable();
		arryScorer.errorTable = adaptor.getErrorsTable();
		arryScorer.ionMap = TheoreticalSpectrumFactory.standardIonMap;
		return arryScorer;
	}
	
	public static ArraySpectrumComparator loadMixtureComparator(String path){
		MixturePeakScoreLearner peakscorer = MixturePeakScoreLearner.loadComparatorLocal(path);
		MixtureScorerAdapter adaptor = new MixtureScorerAdapter(MixTheoSpectrumFactory.mixTypeMap, 
				MixTheoSpectrumFactory.mixIonMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		arryScorer.scoreTable = adaptor.getTable();
		arryScorer.errorTable = adaptor.getErrorsTable();
		arryScorer.ionMap = MixTheoSpectrumFactory.mixIonMap;
		return arryScorer;
	}
	
	
	public static ArraySpectrumComparator loadLinkedComparator(String path){
		LinkedPeptidePeakScoreLearner peakscorer = LinkedPeptidePeakScoreLearner.loadComparator(path);
		LinkedScorerAdapter adaptor = new LinkedScorerAdapter(LinkedTheoSpectrumFactory.linkedTypeMap, 
				LinkedTheoSpectrumFactory.linkedIonMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		arryScorer.scoreTable = adaptor.getTable();
		arryScorer.errorTable = adaptor.getErrorsTable();
		arryScorer.ionMap = LinkedTheoSpectrumFactory.linkedIonMap;
		return arryScorer;
	}
	
	public static void compareScoringMethods(){
		MZXMLReader reader = new MZXMLReader("..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mzXML");
		Spectrum s = reader.getSpectrum(2328);
		s.windowFilterPeaks(8, 25);
		s.computePeakRank();
		Spectrum arry = ArraySpectrum.getRankSpectrum(s);
		Peptide pep = new Peptide("ENEMLAQDK.2");
		TheoreticalSpectrum th = new TheoreticalSpectrum(pep);
		System.out.println("Theoretical spectrum has size: " + th.getPeak().size());
		Map<String, IonType> typeMap = TheoreticalSpectrumFactory.createStandardIonTypeMap();
		IonTypeMapper ionMap = IonTypeMapper.createIonTypeMap(typeMap.values());
		ArrayTheoreticalSpectrum arryTh0 = TheoreticalSpectrumFactory.getArrayTheoSpectrum(pep, typeMap, ionMap);
		Spectrum arryTh = TheoreticalSpectrumFactory.getTheoSpectrumX(pep.getPeptide(), (int)pep.getCharge(), typeMap, ionMap);
		RankBaseScoreLearner peakscorer = RankBaseScoreLearner.loadComparator("../mixture_linked/yeast_single_model_realannotated_win10_25.o");
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		RankBaseScoreAdapter adaptor = new RankBaseScoreAdapter(typeMap, ionMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		arryScorer.scoreTable = adaptor.getTable();
		arryScorer.errorTable = adaptor.getErrorsTable();
		arryScorer.ionMap = ionMap;
		long start = (new GregorianCalendar()).getTimeInMillis();
		double score = 0;
		for(int i = 0; i < 1; i ++){
			score = scorer.compare(th, s);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs score is: " + score);
		start = (new GregorianCalendar()).getTimeInMillis();
		score = 0;
		for(int i = 0; i < 1; i ++){
			score = arryScorer.compare(arryTh, arry);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs score is: " + score);
		
	}
	
	public static void compareLinkedScoringMethods(){
		String file = "..//mixture_linked//lib_sumo1_spectra_with_pyroQ_with_acetylation.mgf";
		SpectrumLib lib = new SpectrumLib(file, "MGF");
		List<Spectrum> specList = lib.getSpectrumList();
		long start = (new GregorianCalendar()).getTimeInMillis();
		Map<String, LinkedIonType> typeMap = LinkedTheoSpectrumFactory.linkedTypeMap;
		LinkedPeptidePeakScoreLearner peakscorer = LinkedPeptidePeakScoreLearner.loadComparator("../mixture_linked/lib_sumo_linked_score_model1.o");
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		LinkedScorerAdapter adaptor = new LinkedScorerAdapter(typeMap, LinkedTheoSpectrumFactory.linkedIonMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		arryScorer.scoreTable = adaptor.getTable();
		arryScorer.errorTable = adaptor.getErrorsTable();
		arryScorer.ionMap = LinkedTheoSpectrumFactory.linkedIonMap;
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			System.out.println("spec is: " + s.peptide + "\t" + s.parentMass + "\t" + s.charge);
			s.windowFilterPeaks(10, 25);
			s.computePeakRank();
			Spectrum arry = ArraySpectrum.getRankSpectrum(s, peakscorer.getRankInterval());
			Peptide pep1 = new Peptide(s.peptide.split("--")[0], 1);
			Peptide pep2 = new Peptide(s.peptide.split("--")[1], 1);
			System.out.println(pep1 + " & " + pep2);
			LinkedPeptide pep = new LinkedPeptide(pep1, pep2, s.charge, pep1.getPeptide().indexOf('K')+1, 6);
			System.out.println("peptide is: " + pep);
			TheoreticalSpectrum th = new TheoreticalSpectrum(pep.getPeptides()[0], pep.getPeptides()[1], (short)s.charge, true);
			System.out.println("Theoretical spectrum has size: " + th.getPeak().size());	
			ArrayTheoreticalSpectrum arryTh = (ArrayTheoreticalSpectrum)LinkedTheoSpectrumFactory.getLinkedTheoSpectrum(pep, typeMap, LinkedTheoSpectrumFactory.linkedIonMap);
			System.out.println("Theoretical array spectrum has size: " + arryTh.getMassIntensityList()[0].length);
			double score = 0, score1 = 0;
			score = scorer.compare(th, s);
			score1 = arryScorer.compare(arryTh, arry);
			System.out.println("matching spectra with scores:\t " + score + "\t" + score1);
			//break;
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000);
	}
	

	public static void compareMixScoringMethods(){
		MZXMLReader reader = new MZXMLReader("..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mzXML");
		Spectrum s = reader.getSpectrum(699);
		s.windowFilterPeaks(8, 25);
		s.computePeakRank();
		ArraySpectrum arry = ArraySpectrum.getRankSpectrum(s);
		String p1 = "VSPADAAK.2";
		String p2 = "WAARKK.2";
		String str1 = "VSPADAAK";
		String str2 = "WAARKK";
		int charge1 = 2;
		int charge2 = 2;
		Peptide pep1 = new Peptide(p1);
		Peptide pep2 = new Peptide(p2);
		TheoreticalSpectrum th = new TheoreticalSpectrum(p1, p2);
		System.out.println("Theoretical spectrum has size: " + th.getPeak().size());
		ArraySpectrum arryTh = (ArraySpectrum)MixTheoSpectrumFactory.getMixTheoSpectrum(str1, str2, charge1, charge2, MixTheoSpectrumFactory.mixTypeMap, MixTheoSpectrumFactory.mixIonMap);
		ArrayTheoreticalSpectrum arryTh1 = (ArrayTheoreticalSpectrum)TheoreticalSpectrumFactory.getTheoSpectrumX(str1, charge1, TheoreticalSpectrumFactory.standardTypeMap, TheoreticalSpectrumFactory.standardIonMap);
		ArrayTheoreticalSpectrum arryTh2 = (ArrayTheoreticalSpectrum)TheoreticalSpectrumFactory.getTheoSpectrumX(str2, charge2, TheoreticalSpectrumFactory.standardTypeMap, TheoreticalSpectrumFactory.standardIonMap);
		//ArraySpectrum combineArry = (ArrayTheoreticalSpectrum)MixTheoSpectrumFactory.getMixTheoSpectrum(arryTh1, arryTh2);
		MixturePeakScoreLearner peakscorer = MixturePeakScoreLearner.loadComparator("../mixture_linked/yeast_simmix_alpha_generic_12_25.o");
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		MixtureScorerAdapter adaptor = new MixtureScorerAdapter(MixTheoSpectrumFactory.mixTypeMap, MixTheoSpectrumFactory.mixIonMap, peakscorer);
		ArraySpectrumComparator arryScorer = new ArraySpectrumComparator(0.5, peakscorer);
		System.out.println("ArrayTheo spect has size: " + arryTh.getMassIntensityList()[0].length);
		arryScorer.scoreTable = adaptor.getTable();
		arryScorer.errorTable = adaptor.getErrorsTable();
		arryScorer.ionMap = MixTheoSpectrumFactory.mixIonMap;
		long start = (new GregorianCalendar()).getTimeInMillis();
		double score = 0;
		for(int i = 0; i < 1; i ++){
			score = scorer.compare(th, s);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs score is: " + score);
		start = (new GregorianCalendar()).getTimeInMillis();
		score = 0;
		for(int i = 0; i < 1; i ++){
			score = arryScorer.compare(arryTh, arry);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs score is: " + score);
		
		for(int i = 0; i < 1000000; i ++){
			ArraySpectrum combineArry = (ArrayTheoreticalSpectrum)MixTheoSpectrumFactory.getMixTheoSpectrum(arryTh1, arryTh2);
			score = arryScorer.compare(combineArry, arry);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs score is: " + score);

	}
	
	public static void main(String[] args){
		//compareScoringMethods();
		//compareMixScoringMethods();
		compareLinkedScoringMethods();	
	}
}
