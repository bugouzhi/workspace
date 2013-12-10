package org.Spectrums;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * create a map for spectrumlib so spectrum in the lib
 * can be lookup using some property like scan number
 * @author Jian
 *
 */
public class SpectrumLibMap extends SpectrumLib{
	private Map<Integer, Spectrum> scanMap;
	
	public SpectrumLibMap(String file, String format){
		super(file, format);
		createScanMap();
	}
	public void createScanMap(){
		this.scanMap = new HashMap<Integer, Spectrum>();
		for(Iterator it = this.getAllSpectrums().iterator(); it.hasNext();){
			Spectrum s = (Spectrum)it.next();
			this.scanMap.put(s.scanNumber, s);
		}
	}
	
	public Spectrum getSpecByScan(int scan){
		return this.scanMap.get(scan);
	}
}
