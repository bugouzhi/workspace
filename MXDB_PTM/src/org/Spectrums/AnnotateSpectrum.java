package org.Spectrums;
/**\
 * Annotate a spectrum
 * @author Jian Wang
 *
 */
public class AnnotateSpectrum {
	public  static void annotateSpectrum(){
		String spectrumFile ="../mixture_linked//msdata/philAndrews/zymogen/zgm_xlink_BS3_032803_scx10.mzXML";
		int scan = 3471;
		String peptide = "SPPNENIIINNPSRPWWER";
		int charge = 2;
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		Spectrum s = reader.getSpectrum(scan);
		s.computePeakRank();
		Peptide p = new Peptide(peptide, charge);
		System.out.println("peptide: " + p + "\t" + p.getParentmass());
		System.out.println("spectrum has charge: " + s.charge);
		TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		TheoreticalSpectrum t = new TheoreticalSpectrum(p);
		t.analyzeAnnotation(s, peptide);
	}
	
	public static void main(String[] args){
		annotateSpectrum();
	}
}
