package org.Spectrums;

import java.util.List;

/**
 * Generate theoretical spectrum for fragment ions rather than normal ions
 * @author Jian Wang
 *
 */
public class TheoreticalFragmentSpectrum {
	/**
	 * Generate theoretical fragment masses for a fragments
	 */
	
	public static TheoreticalSpectrum getTheoreticalSpectrum(Peptide p, LabelledPeak matchedPeak){
		int offSet = matchedPeak.getPos();
		String pep = "";
		Peptide sub = null;
		TheoreticalSpectrum t = null;
		if(matchedPeak.isPrefixPeak()){
			pep = sub.getPeptide();
			sub =  new Peptide(pep.substring(0, offSet));		
		}else if(matchedPeak.isSuffixPeak()){
			pep = p.getPeptide();
			sub = new Peptide(pep.substring(pep.length()-offSet, pep.length()));
		}
		Peptide subPep = new Peptide(sub);
		return getSubPeptideSpectrum(subPep, matchedPeak);
	}
	
	private static TheoreticalSpectrum getSubPeptideSpectrum(Peptide sub, LabelledPeak matchedPeak){
		TheoreticalSpectrum t = new TheoreticalSpectrum(sub);
		return getSubPeptideSpectrum(t, matchedPeak);
		
	}
	
	private static TheoreticalSpectrum getLinkedSubPeptideSpectrum(Peptide sub, LabelledPeak matchedPeak, int linkedCharge){
		TheoreticalSpectrum t = new TheoreticalSpectrum(sub, linkedCharge);
		return getSubPeptideSpectrum(t, matchedPeak);
	}
	
	/**
	 * modify the masses of prefix/suffix peak accordingly
	 * @param t
	 * @return
	 */
	private static TheoreticalSpectrum getSubPeptideSpectrum(TheoreticalSpectrum t, LabelledPeak matchedPeak){
		if(matchedPeak.isPrefixPeak()){
			List<Peak> pList = t.getPeaks();
			for(int i = 0; i < pList.size(); i++){
				LabelledPeak current = (LabelledPeak)pList.get(i);
				if(current.isSuffixPeak()){
					current.setMoz(current.getMass()
							+(Mass.getIonMod(matchedPeak.getPeakType())
							-Mass.PROTON_MASS)/matchedPeak.getCharge());
				}
			}
		}
		
		if(matchedPeak.isSuffixPeak()){
			List<Peak> pList = t.getPeaks();
			for(int i = 0; i < pList.size(); i++){
				LabelledPeak current = (LabelledPeak)pList.get(i);
				if(current.isPrefixPeak()){
					current.setMoz(current.getMass()
							+(Mass.getIonMod(matchedPeak.getPeakType())
							-Mass.PROTON_MASS)/matchedPeak.getCharge());
				}
			}
			
		}
		return t;
	}
	
	
	public static TheoreticalSpectrum getTheoreticalSpectrum(LinkedPeptide p, LabelledPeak matchedPeak){
		int offSet = matchedPeak.getPos();
		String pep = "";
		Peptide sub = null;
		Peptide full = null;
		TheoreticalSpectrum t = null;
		Peptide[] peptides = p.peptides;
		if(matchedPeak.getPep().getPeptide().equals(peptides[0].getPeptide())){
			sub = peptides[0];
			full = peptides[1];
		}
		if(matchedPeak.getPep().getPeptide().equals(peptides[1].getPeptide())){
			sub = peptides[1];
			full = peptides[0];
		}
		if(matchedPeak.isPrefixPeak()){
			pep = sub.getPeptide();
			sub = new Peptide(sub);
			sub.setPeptide(pep.substring(0, offSet));
		}else if(matchedPeak.isSuffixPeak()){
			pep = sub.getPeptide();
			sub = new Peptide(sub);
			sub.setPeptide(pep.substring(pep.length()-offSet, pep.length()));
		}
		full.setPtmmasses(new double[]{sub.getParentmass()});
		TheoreticalSpectrum t1 = getLinkedSubPeptideSpectrum(sub, matchedPeak, p.getCharge());
		TheoreticalSpectrum t2 = new TheoreticalSpectrum(full, matchedPeak.getPep().getCharge());
		TheoreticalSpectrum linked = new TheoreticalSpectrum(full, sub, p.getCharge(), false);
		return linked;
	}
	
	
	
	
}
