package org.Spectrums;


import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import org.systemsbiology.jrap.stax.*;
public class SortedMZXMLReader extends MZXMLReader{
	private int currentIndex = 0;
	private Spectrum current;
	private List<Integer> scanList;
	private int minScan = 0;
	private int maxScan = Integer.MAX_VALUE;
	
	public SortedMZXMLReader(String filename){
		super(filename);
		initialize();
	}
	
	//only read spectrum within certain range
	public SortedMZXMLReader(String filename, int minScan, int maxScan){
		super(filename);
		this.minScan = minScan;
		this.maxScan = maxScan;
		initialize();
	}
	
	public void initialize(){
		TreeMap<Double, Integer> scanMap = new TreeMap<Double, Integer>();
		for(int i = this.minScan; i < this.parser.getScanCount() && i < this.maxScan; i++){
			Scan s = parser.rap(i);
			if(s!= null && s.getHeader().getMsLevel() == 2){
				scanMap.put(new Double(s.getHeader().getPrecursorMz()+ Math.random()* 0.000001), s.getHeader().getNum());
			}
		}
		this.scanList = new ArrayList<Integer>(scanMap.keySet().size());
		for(Iterator<Double> iter = scanMap.navigableKeySet().iterator(); iter.hasNext();){
			this.scanList.add(scanMap.get(iter.next()));
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
			return this.getSpectrum(scanList.get(currentIndex++));
		}else{
			return null;
		}
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}

}
