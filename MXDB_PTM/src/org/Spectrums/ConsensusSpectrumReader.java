package org.Spectrums;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
public class ConsensusSpectrumReader implements Iterator<Spectrum>{
	private MZXMLReader reader;
	private List<Spectrum> specList;
	public int numNeighbors=3;
	private Spectrum  current;
	public double minInt=0;
	public boolean decoyMode = false;
	
	public ConsensusSpectrumReader(String spectrumFile){
		this.reader = new MZXMLReader(spectrumFile);
		this.specList = new ArrayList<Spectrum>();
	}
	
	public ConsensusSpectrumReader(MZXMLReader reader){
		this.reader = reader;
		this.specList = new ArrayList<Spectrum>();
	}
	
	public Spectrum getSpectrum(int scan){
		Spectrum currentSpec = this.reader.getSpectrum(scan);
		getNeighborScans(currentSpec);
		currentSpec.filterPeaksByIntensity(this.minInt);
		//currentSpec.filterPeaks(1000);
		computeConsensus(currentSpec, this.specList);
		return currentSpec;
	}
	
	private Spectrum getConsensusSpectrum(){
		this.current = this.reader.next();
		getNeighborScans(this.current);
		//this.current.filterPeaksByIntensity(this.minInt);
		this.current.filterPeaks(1000);
		computeConsensus(this.current, this.specList);
		return current;
		
	}
	
	private void getNeighborScans(Spectrum s){
		int minScan = 0;
		int maxScan = this.reader.getSpectrumCount();
		int Cycle = 35;
		if(decoyMode)  //for decoy mode we pick a random scan rather than correct +/- SWATH scans
			Cycle = 10;//(int)Math.floor(Math.random()*25)+5;
		int currentScan = s.scanNumber;
		//System.out.println("current Scan: " + currentScan);
		this.specList.clear();
		for(int i = 1; i <= this.numNeighbors; i++){
			int scan = -1*i*Cycle + currentScan;
			if(scan > minScan && scan < maxScan){
				Spectrum neigh = this.reader.getSpectrum(scan);
				//System.out.println("gettting scan: " + scan + "\tmz: " + neigh.parentMass);
				this.specList.add(neigh);
			}
		}
		
		for(int i = 1; i <= this.numNeighbors; i++){
			int scan = 1*i*Cycle + currentScan;
			if(scan > minScan && scan < maxScan){
				Spectrum neigh = this.reader.getSpectrum(scan);
				//System.out.println("getting scan: " + scan + "\tmz: " + neigh.parentMass);
				this.specList.add(neigh);
			}
		}
	}
	
	private void computeConsensus(Spectrum s, List<Spectrum> specList){
		for(int i = 0; i < specList.size(); i++){
			Spectrum current = specList.get(i);
			//current.filterPeaksByIntensity(this.minInt);
			current.filterPeaks(1000);
			current.computeConsensus(s, 0.05);
		}
	}
	
	
	
	@Override
	public boolean hasNext() {
		return this.reader.hasNext();
	}

	@Override
	public Spectrum next() {
		return getConsensusSpectrum();
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}
	
	public static  void testConsensus(){
		String spectrumFile = "../mixture_linked/msdata/gringar/swath_development/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		ConsensusSpectrumReader reader = new ConsensusSpectrumReader(spectrumFile);
		reader.numNeighbors = 2;
		reader.minInt = 0;
		reader.decoyMode = true;
		while(reader.hasNext()){
			Spectrum s = reader.next();
			//s = reader.getSpectrum(11624);
			if(s.scanNumber != 11624){
				//continue;
			}
			//System.out.println("consensus: " + s.scanNumber);
			//s.filterPeaksByIntensity(10);
			//s.windowFilterPeaks2(15, 25);
			s.filterPeaksByRankScore(4);
			//s.windowFilterPeaks2(15, 25);
			//System.out.println(s);
			//return;
		}
	}
	
	public static void main(String[] args){
		testConsensus();
	}

}
