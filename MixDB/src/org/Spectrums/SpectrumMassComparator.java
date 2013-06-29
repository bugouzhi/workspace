package org.Spectrums;
import java.util.Comparator;

public class SpectrumMassComparator implements Comparator{
	public static SpectrumMassComparator comparator = new SpectrumMassComparator();
	
	@Override
	public int compare(Object arg0, Object arg1) {
		Spectrum s1 = (Spectrum)arg0;
		Spectrum s2 = (Spectrum)arg1;
		if(s1.parentMass > s2.parentMass){
			return 1;
		}else if(s1 == s2){
			return 0;
		}else{
			return -1;
		}
	}

}
