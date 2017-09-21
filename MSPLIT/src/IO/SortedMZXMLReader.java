package IO;
/**
 * Reading mzxml file in order according to the precursor mass of the spectrum
 */


import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import org.systemsbiology.jrap.stax.*;

import Spectrum.Spectrum;
public class SortedMZXMLReader extends MZXMLReader{
	private int currentIndex = 0;
	private Spectrum current;
	private List<Integer> scanList;
	private List<Integer> indexList;
	
	public SortedMZXMLReader(String filename){
		super(filename);
		initialize();
	}
	
	
	public void initialize(){
		TreeMap<Double, Integer> scanMap = new TreeMap<Double, Integer>();
		TreeMap<Double, Integer> indexMap = new TreeMap<Double, Integer>();
		int index = 0;
		for(int i = 1; i < this.getParser().getMaxScanNumber(); i++){
			Scan s = parser.rap(i);
			if(s!= null && s.getHeader().getMsLevel() == 2){
				index++;
				Double key = s.getHeader().getPrecursorMz()+ Math.random()* 0.000001;
				scanMap.put(key, s.getHeader().getNum());
				indexMap.put(key, index);
			}
		}
		this.scanList = new ArrayList<Integer>(scanMap.keySet().size());
		this.indexList = new ArrayList<Integer>(scanMap.keySet().size());
		for(Iterator<Double> iter = scanMap.navigableKeySet().iterator(); iter.hasNext();){
			Double key = iter.next();
			this.scanList.add(scanMap.get(key));
			this.indexList.add(indexMap.get(key));
		}
		System.out.println("read in spectra from mzXML: " + this.scanList.size());
	}
	
	@Override
	public boolean hasNext() {
		return this.scanList != null && currentIndex < this.scanList.size();
	}

	@Override
	public Spectrum next() {
		if(hasNext()){
			Spectrum s = this.getSpectrum(scanList.get(currentIndex));
			s.specIndex = this.indexList.get(currentIndex++);
			return s;
		}else{
			return null;
		}
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}

}
