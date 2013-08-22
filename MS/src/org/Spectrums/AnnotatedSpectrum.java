package org.Spectrums;

import java.util.HashMap;
import java.util.Map;

/**\
 * A spectrum allowing for various annotation 
 * The annotation are stored as a general map of the annotation name and their values
 * @author Jian Wang
 *
 */
public class AnnotatedSpectrum extends Spectrum{
	private Map<String, Object> annotation;
	
	public AnnotatedSpectrum(){
		super();
		this.annotation = new HashMap<String, Object>();
	}
	
	public Map<String, Object> getAnnotation() {
		return annotation;
	}

	public void setAnnotation(Map<String, Object> annotation) {
		this.annotation = annotation;
	}

	public  static void annotateSpectrum(){
//		String spectrumFile ="../mixture_linked//msdata/philAndrews/zymogen/zgm_xlink_BS3_032803_scx10.mzXML";
//		int scan = 3471;
//		String peptide = "SPPNENIIINNPSRPWWER";
//		int charge = 2;
//		MZXMLReader reader = new MZXMLReader(spectrumFile);
//		Spectrum s = reader.getSpectrum(scan);
//		s.computePeakRank();
//		Peptide p = new Peptide(peptide, charge);
//		System.out.println("peptide: " + p + "\t" + p.getParentmass());
//		System.out.println("spectrum has charge: " + s.charge);
//		TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
//		TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
//		TheoreticalSpectrum t = new TheoreticalSpectrum(p);
//		t.analyzeAnnotation(s, peptide);
	}
	
	public static void main(String[] args){
		annotateSpectrum();
	}
}
