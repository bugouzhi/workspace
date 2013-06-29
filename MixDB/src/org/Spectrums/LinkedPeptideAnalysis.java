package org.Spectrums;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import Utils.FileIOUtils;

/**
 * Some utility function for analyzing linked peptides and their annotation
 * @author Jian Wang
 *
 */
public class LinkedPeptideAnalysis {
	public static void getIonsRank(Spectrum s){
		String[] peptides = s.peptide.split(" & ");
		s.computePeakRank();
	}
	
	public static void analyzeSumoSpectra(){
		String file = "..\\mixture_linked\\lib_sumo2_spectra_crosslink_search_lib_sumo_linked_model_win8_25_floatK.mgf";
		SpectrumLib lib = new SpectrumLib(file, "MGF");
		for(Iterator<Spectrum> it = lib.getAllSpectrums().iterator(); it.hasNext(); ){
			Spectrum current = it.next();
			current.windowFilterPeaks(8, 25);
			System.out.println("after filtering we have peaks " + current.getPeak().size());
			current.removePrecursors(0.5);
			current.computePeakRank();
			System.out.println(current.spectrumName + "\t" + current.charge + "\t" + current.parentMass 
					+ "\tfiltered spectrum has peaks:\t" + current.getPeak().size());
			String[] peptides = current.peptide.split("--");
			LinkedPeptide lp = new LinkedPeptide(current.peptide, current.charge, 5, 6);
			TheoreticalSpectrum linkedSpect = new TheoreticalSpectrum(lp.peptides[0], lp.peptides[1], lp.getCharge(), false);
			double[] stat = linkedSpect.analyzeMixtureAnnotation(current, peptides[0], peptides[1], 0.5);
			System.out.println(current.spectrumName + "\t" + current.getPeptide() + "\t" +  current.parentMass + "\tbest\t" +  lp + "\t" + lp.getParentmass() + "\t" + lp.getCharge() + "\t" 
					+ " with score:\t" +  "\t" + stat[0] + "\t" + stat[1] + "\t"
					+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
					+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
					+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12]);
		}
	}
	
	public static void analyzeSumoAnnotation(String spectrumFile, String annotationFile){
		//String spectrumLibFile = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mzXML";
		List<String> lines = FileIOUtils.createListFromFile(annotationFile);
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		for(int i = 0; i < lines.size(); i++){
			String[] tokens = lines.get(i).split("\\s+");
			int scan = Integer.parseInt(tokens[2]);
			Spectrum s = reader.getSpectrum(scan);
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			//int charge = Integer.parseInt(tokens[5]);
			int charge = s.charge;
			String peptide = tokens[3];
			System.out.println(s.spectrumName + "\t peptide is: " + peptide + " charge is : " + s.charge);
			LinkedPeptide lp = new LinkedPeptide(peptide, charge, 1, 6);
			TheoreticalSpectrum linkedSpect = new TheoreticalSpectrum(lp, lp.getCharge(), false);
			SpectrumComparator scorer2 = SpectrumUtil.getLinkedPeptideScorer("..\\mixture_linked\\testScoring.o");
			MixtureSpectrumScorer.detail = true;
			scorer2.compare(linkedSpect, s);
//			double[] stat = linkedSpect.analyzeMixtureAnnotation(s, peptides[0], peptides[1], 0.5);
//			System.out.println("Spectrum: " + s.getPeptide() + "\t" +  s.parentMass + "\tbest\t" +  lp + "\t" + lp.getParentmass() + "\t" + lp.getCharge() + "\t" 
//					+ " with score:\t" +  "\t" + stat[0] + "\t" + stat[1] + "\t"
//					+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
//					+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
//					+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12]);
		}
	}

	public static void getLinkedModel(){
		//SpectrumComparator scorer1 = SpectrumUtil.getRankBaseScorer("..\\MSPLib\\Lib\\yeast.msp");
		TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		MixtureSpectrumScorer scorer2 = (MixtureSpectrumScorer)SpectrumUtil.getLinkedPeptideScorer("..\\mixture_linked\\lib_sumo_spectra_training_part1.mgf");
		((LinkedPeptidePeakScoreLearner)scorer2.comp).writeLibToFile("..\\mixture_linked\\test_model.o");
	}
	
	public static void getLinkedSinglePeptideModel(){
		TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		SimpleProbabilisticScorer scorer2 = (SimpleProbabilisticScorer)SpectrumUtil.getLPeakRankBaseScorer("..\\MSPLib\\Lib\\yeast.msp");
		((LPeakRankBaseScorer)scorer2.comp).writeLibToFile("..\\mixture_linked\\yeast_single_peptide_score_model_unfilter.o");
	}
	
	public static void getSinglePeptideModel(){
		//TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		//TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		//SimpleProbabilisticScorer scorer2 = (SimpleProbabilisticScorer)SpectrumUtil.getRankBaseScorer("..\\MSPLib\\Lib\\yeast.msp");
		SimpleProbabilisticScorer scorer2 = (SimpleProbabilisticScorer)SpectrumUtil.getRankBaseScorer("..\\mixture_linked\\human_heck_inspecttraining.mgf");
		((RankBaseScoreLearner)scorer2.comp).writeLibToFile("..\\mixture_linked\\test.o");
	}
	
	public static void getLinkedMixtureModel(){
		//TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		//TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		MixtureSpectrumScorer scorer2 = (MixtureSpectrumScorer)SpectrumUtil.getLMixtureScorer("..\\mixture_linked\\mixtures100000_alpha0.3.mgf");
		((LMixturePeakScoreLearner)scorer2.comp).writeLibToFile("..\\mixture_linked\\mixtures_alpha0.3_models.o");
	}
	
	public static void getMixtureModel(){
		//TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		//TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		MixtureSpectrumScorer scorer2 = (MixtureSpectrumScorer)SpectrumUtil.getMixtureScorer("..\\mixture_linked\\mixtures100000_generic.mgf");
		((MixturePeakScoreLearner)scorer2.comp).writeLibToFile("..\\mixture_linked\\yeast_simmix_alpha_generic_8_25.o");
	}
	public static void main(String[] args){
		//String spectrumFile = "..\\mixture_linked\\Linked_peptide_library\\sumo_lib\\20101008_Sumo_Library_4349_Bo.mzXML";
		//String annotationFile = "..\\mixture_linked\\test0.txt";
		//analyzeSumoSpectra();
		//getLinkedModel();
		//getLinkedSinglePeptideModel();
		//getLinkedMixtureModel();
		getMixtureModel();
		//analyzeSumoAnnotation(spectrumFile, annotationFile);
		//getSinglePeptideModel();
	}
}
