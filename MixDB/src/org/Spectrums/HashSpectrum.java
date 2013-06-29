package org.Spectrums;
import java.util.Map;
import org.apache.commons.collections.MultiMap;
import org.apache.commons.collections.map.MultiValueMap;

public class HashSpectrum extends Spectrum{
	Map peakTable;
	public HashSpectrum(Spectrum s){
		this.charge = s.charge;
		this.spectrumName = s.spectrumName;
		this.parentMass = s.parentMass;
		this.peakTable = new MultiValueMap();
		for(int i = 0; i < s.getPeak().size(); i++){
			Peak current = s.getPeak().get(i);
			this.peakTable.put(Math.round(current.getMass()), current);
		}
	}

}
