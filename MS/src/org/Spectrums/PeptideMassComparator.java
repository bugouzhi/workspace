package org.Spectrums;
import java.util.Comparator;
/**
 * Impose ordering of peptides according to their masses
 * note that in order to allow the notion that different peptides
 * can has exactly the same masses, when mass is equal the comparator do 
 * not return zero, which is reserve for cases when the two peptides are the same
 * @author Jian Wang
 *
 */
public class PeptideMassComparator implements Comparator{
public static PeptideMassComparator comparator = new PeptideMassComparator();
	@Override
	public int compare(Object arg0, Object arg1) {
		Peptide s1 = (Peptide)arg0;
		Peptide s2 = (Peptide)arg1;
		if(s1.getParentmass() >= s2.getParentmass()){
			return 1;
		}else if(s1 == s2){
			return 0;
		}else{
			return -1;
		}
	}
}
