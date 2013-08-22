package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

/**
 * Contain various method that analyze mixtrure spectrum
 * @author jian wang
 *
 */
public class MixturePeptideAnalysis {
	public static void analyzeMixtureIds(){
		String inputfile = "..\\mixture_linked\\Linked_peptide_library\\sumo_lib\\20101008_Sumo_Library_4351_Bo.mzXML";
		MZXMLReader iterator = new MZXMLReader(inputfile);
		Spectrum query = iterator.getSpectrum(4119);
		query.windowFilterPeaks(10, 25);
		query.computePeakRank();
		String peptide1 = "QQQTGGAWKMETPFRAK.2";
		String peptide2 = "QQQTGGAYKHETDFRAK.2";
		TheoreticalSpectrum t = new TheoreticalSpectrum(peptide1, peptide2);
		double[] stat =t.analyzeMixtureAnnotation(query, peptide1, peptide2);
		double score1 = 0;
		double score2 = 0;
		SpectrumComparator scorer1 = SpectrumUtil.getRankBaseScorer("..\\MSPLib\\Lib\\Ecoli.msp");
		List<Spectrum> list = new ArrayList();
		list.add(new TheoreticalSpectrum(peptide1));
		list.add(new TheoreticalSpectrum(peptide2));
		list.add(new TheoreticalSpectrum("EIIDGFLKFQREAFPK.2"));		
		SpectrumLibSearcher searcher = new SpectrumLibSearcher(list, scorer1, scorer1);
		searcher.bestCandidates(query, 10);
//		System.out.print("Spectrum: " + query.getPeptide() + " has best match: " +  peptide1 + " & " + peptide2
//				+ " with score:\t" +  score1 + "\t" + score2 + "\t"
//				+ "\t" + stat[0] + "\t" + stat[1] + "\t"
//				+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
//				+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
//				+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12]);
//		//System.out.println("\t" + checkPeptidepair(bestpeptide, query.peptide));	
		System.out.println();
	}
	
	public static void main(String args[]){
		analyzeMixtureIds();
	}
}
