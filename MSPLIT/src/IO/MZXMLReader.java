package IO;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Iterator;
import java.util.TreeMap;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreeNode;

import Spectrum.Mass;
import Spectrum.Peak;
import Spectrum.Spectrum;
//import org.Spectrums.SpectrumUtil;
import org.systemsbiology.jrap.stax.*;

/**
 * Reader of mzxml files
 * @author Jian
 *
 */
//a reader that can read mzXML files
public class MZXMLReader implements SpectrumReader{
	protected String filename;
	public MSXMLParser parser;
	private int currentScan = 1;
	private Spectrum current;
	private int maxScan = 0;
	public MZXMLReader(String filename){
		this.filename = filename;
		this.parser = new MSXMLParser(this.filename);
		this.maxScan = this.parser.getMaxScanNumber();
		//System.out.println("Total Number of spectra in file: " + this.parser.getScanCount());
		//System.out.println("Beginning scan " + currentScan);
		next();
	}
	
	public MZXMLReader(String filename,int startScan){
		this.filename = filename;
		this.parser = new MSXMLParser(this.filename);
		this.maxScan = this.parser.getMaxScanNumber();
		this.currentScan = startScan-1;
		//System.out.println("Total Number of spectra in file: " + this.parser.getScanCount());
		//System.out.println("Beginning scan " + currentScan);
		next();
	}
	
	public MZXMLReader(MSXMLParser parser){
		this.parser = parser;
	}
	
	
	public List<Spectrum> readAllMS2Spectra(){
		List<Spectrum> specList = new ArrayList<Spectrum>();
		for(int i = 1; i < this.parser.getScanCount(); i++){
			List<Peak> peakList = new ArrayList<Peak>();
			Scan current = parser.rap(i);
			if(current.getHeader().getMsLevel() >= 2){
				Spectrum s = getSpectrum(i);
				if(s.charge < 0){ //fix no-charge information artifect temporarily
					s.charge = 0;
				}
				specList.add(s);
			}			
		}
		return specList; 
	}
	
	public MSXMLParser getParser() {
		return parser;
	}
	
	public int getSpectrumCount(){
		return this.parser.getMaxScanNumber();
	}
	
	public TreeModel getStructuredSpectrum(){
		DefaultMutableTreeNode node = new DefaultMutableTreeNode(this.filename); //use name of file as root object
		TreeModel tree = new DefaultTreeModel(node);
		int currentLevel = 0;
		DefaultMutableTreeNode currentParent=node;
		for(int i = 1; i < this.parser.getScanCount(); i++){
			Scan current = parser.rap(i);
			Spectrum s = getSpectrum(current);
			DefaultMutableTreeNode newNode = new DefaultMutableTreeNode(s);
			if(current.getHeader().getMsLevel() > currentLevel){				
				currentParent.add(newNode);
			}else if(current.getHeader().getMsLevel() == currentLevel){
				currentParent = (DefaultMutableTreeNode)currentParent.getParent();
				currentParent.add(newNode);
				
			}else{
				int j = currentLevel;
				while(j >= current.getHeader().getMsLevel()){
					currentParent = (DefaultMutableTreeNode)currentParent.getParent();
					j--;
				}
				currentParent.add(newNode);
			}
			currentLevel = current.getHeader().getMsLevel();
			currentParent = newNode;
		}
		return tree;
	}
	
	public int[] getSpectrumStat(){
		int count1 = 0, count2 = 0, count3 = 0;
		for(int i = 1; i < this.parser.getScanCount(); i++){
			Scan current = parser.rap(i);
			if(current.getHeader().getMsLevel() == 1){
				count1++;
			}
			if(current.getHeader().getMsLevel() == 2){
				count2++;
			}
			if(current.getHeader().getMsLevel() >= 3){
				count3++;
			}
		}
		//System.out.println("Total MS1: " + count1 + "\tMS2: " + count2 + "\tMS3: " + count3);
		return new int[]{count1, count2, count3};
	}
	
	public void setParser(MSXMLParser parser) {
		this.parser = parser;
	}

	public Spectrum getSpectrum(int scanNum){
		List<Peak> peakList = new ArrayList<Peak>();
		Scan current = this.parser.rap(scanNum);
		if(current == null){
			return null;
		}
		//double pif = this.getPIF(scanNum);
		//System.out.println("Scan: "  + scanNum + "\t" + current.getHeader().getPrecursorMz() + "\t" + current.getHeader().getPrecursorCharge() + "\tpif:\t" + pif);
		Spectrum s = new Spectrum();
		s.parentMass  = current.getHeader().getPrecursorMz();
		s.charge = current.getHeader().getPrecursorCharge();
		if(s.charge < 0 && current.getHeader().getMsLevel() > 1){ //cannot obtain charge information from file we try to compute them from MS1
		//	s.charge = -1*this.getPrecursorCharge(current); 
		}
		double[][] peaks = current.getMassIntensityList();
		for(int j = 0; j < peaks[0].length; j++){
			Peak p = new Peak(peaks[0][j], peaks[1][j]);
			peakList.add(p);
		}
		s.setPeaks(peakList);
		s.spectrumName = "Scan Number: " + scanNum + " Retention Time: " + current.getHeader().getRetentionTime();
		s.scanNumber = scanNum;
		s.rt = getRT(current.getHeader().getRetentionTime());
		//s.upperBound = pif;
		//s.upperBound = current.getHeader().getTotIonCurrent();
		//s.score = this.getPIF(s.scanNumber, s.parentMass, 50, -2.0, 2.0, 1, s.charge)[0];
		return s;
	}
		
	public Map<Integer, Double> getRTScanMapping(){
		Map<Integer,Double> rtMapping = new HashMap<Integer, Double>();
		for(int i = 1; i < this.parser.getScanCount(); i++){
			Scan current = parser.rap(i);
			double rt = getRT(current.getHeader().getRetentionTime());
			int scan = current.getHeader().getNum();
			rtMapping.put(scan, rt);
		}
		return rtMapping;
	}
	
	public Map<Double, Integer> getRTScanMappingReverse(){
		Map<Double, Integer> rtMapping = new TreeMap<Double, Integer>();
		for(int i = 1; i < this.parser.getScanCount(); i++){
			Scan current = parser.rap(i);
			if(current.getHeader().getMsLevel() > 0){ //only map MS1 scan in this direction
				double rt = getRT(current.getHeader().getRetentionTime());
				int scan = current.getHeader().getNum();
				rtMapping.put(rt, scan);
			}
		}
		return rtMapping;
	}
	
	public Map<Double, Spectrum> getSpectrumTimeMap(){
		Map<Double, Spectrum> rtMapping = new TreeMap<Double, Spectrum>();
		for(int i = 1; i < this.parser.getScanCount(); i++){
			Spectrum s = getSpectrum(i);
			rtMapping.put(s.rt, s);
		}
		return rtMapping;
	}

	private static double getRT(String rt){
		rt = rt.replaceAll("[^0-9,\\.]", "");
		return Double.parseDouble(rt);
	}
	
	public static Spectrum getSpectrum(Scan current){
		List<Peak> peakList = new ArrayList<Peak>();
		Spectrum s = new Spectrum();
		s.parentMass  = current.getHeader().getPrecursorMz();
		s.charge = current.getHeader().getPrecursorCharge();
		s.rt = getRT(current.getHeader().getRetentionTime());
//		if(s.charge < 0){ //cannot obtain charge information from file we try to compute them from MS1
//			s.charge = this.getPrecursorCharge(current); 
//		}
		double[][] peaks = current.getMassIntensityList();
		for(int j = 0; j < peaks[0].length; j++){
			Peak p = new Peak(peaks[0][j], peaks[1][j]);
			peakList.add(p);
		}
		s.setPeaks(peakList);
		s.spectrumName = "Scan Number: " + current.getHeader().getNum();
		s.scanNumber = current.getHeader().getNum();
		return s;
	}
	
	public int getPrecursorCharge(Scan spectrum){
		return getPrecursorCharge(spectrum, 70);
	}
	
	public static double massDiff(double mass1, double mass2, int mode){
		if(mode == 1){
			return Math.abs(mass1 - mass2);
		}else{
			return Math.abs((mass1 - mass2)*1000000 / mass2); 
		}

	}
	public int getPrecursorCharge(Scan spectrum, double tolerance){
		//System.out.println("processing MS2: " + spectrum.getNum());
		//System.out.println("MS1 is : " + getPrevScan1(spectrum.getNum()));
		Scan MS1 = this.parser.rap(getPrevScan1(spectrum.getHeader().getNum()));
		double precursorMass = spectrum.getHeader().getPrecursorMz();
		//System.out.println("recorded mass is: " + precursorMass);
		double precursorIntensity = spectrum.getHeader().getPrecursorIntensity();
		double[][] peaks = MS1.getMassIntensityList();
		int precursorIndex = findPrecursor(peaks, precursorMass);	
		double diff = Math.abs(Mass.massDiff(peaks[0][precursorIndex],  precursorMass, 2));
		if(diff < tolerance){
			System.out.println("Scan " + spectrum.getHeader().getNum() + ": matched ms1 peak is: " + peaks[0][precursorIndex] + "\t" + peaks[1][precursorIndex]);
			double[][] expectedProfile = new double[5][3];
			double firstmass = (double)peaks[0][precursorIndex];
			for(int i = 0; i < expectedProfile.length; i++){
				for(int k = 0; k < expectedProfile[0].length; k++){
					expectedProfile[i][k] = firstmass + (1.00335/(i+1))*k;
					//System.out.println("expected: " + expectedProfile[i][k]);
				}
			}
			double[] matchCount = new double[expectedProfile.length];
			for(int i = 0; i < expectedProfile.length; i++){
				int k=1;
				for(int j = precursorIndex+1; j < peaks[0].length-1; j++){			
					if(Mass.massDiff(expectedProfile[i][k], peaks[0][j],2) < tolerance){
						matchCount[i]+= peaks[1][j];
						k++;
					}
					
					if(k >= expectedProfile[i].length){
						break;
					}
				}
			}
			
			int maxInd = -1;
			double max = 0;
			System.out.print("Scan " + spectrum.getHeader().getNum() + "\t: matched profile:\t");
			for(int i = 0; i < matchCount.length; i++){
				System.out.print(matchCount[i] + "\t");
				maxInd = max >= matchCount[i] ? maxInd : i;
				max = max >= matchCount[i] ? max : matchCount[i];
			}
			System.out.println();
			return maxInd+1;
		}	
		return 0;
	}
	public int getPrecursorCharge(int scanNum){
		Scan MS2 = this.parser.rap(scanNum);
		return getPrecursorCharge(MS2);
	}
	
	public int findPrecursor(double[][] peaks, double precursorMass){
		double min = Float.MAX_VALUE;
		int minIndex = -1;
		for(int i = 0; i < peaks[0].length; i++){
				double diff = Math.abs(Mass.massDiff(peaks[0][i],  precursorMass, 2));
				minIndex = min < diff ? minIndex : i; 
				min = min < diff ? min : diff;
		}
		return minIndex;
	}
	
	public boolean isIsotopeCoded(int scanNum, double massShift, double tolerance){
		Scan spectrum = this.parser.rap(scanNum);
		Scan MS1 = this.parser.rap(getPrevScan1(spectrum.getHeader().getNum()));
		double precursorMass = spectrum.getHeader().getPrecursorMz();
		double precursorIntensity = spectrum.getHeader().getPrecursorIntensity();
		double[][] peaks = MS1.getMassIntensityList();
		int precursorIndex = findPrecursor(peaks, precursorMass);	
		double diff = Math.abs(Mass.massDiff(peaks[0][precursorIndex],  precursorMass, 2));
		if(diff < 10){
					System.out.println("Scan number: " + scanNum + " matched ms1 peak is: " + peaks[0][precursorIndex] + "\t" + peaks[1][precursorIndex] + "\t" + diff);
					for(int j = precursorIndex+1; j < peaks[0].length-1; j++){
						if(Mass.massDiff(peaks[0][j],  precursorMass + massShift, 2) < 10){
							System.out.println("Scan Number " + scanNum + " has charge 1 isotope");
							return true;
						}
						if(Mass.massDiff(peaks[0][j],  precursorMass + massShift/2, 2) < 10){
							System.out.println("Scan Number " + scanNum + " has charge 2 isotope");
							return true;
						}
						if(Mass.massDiff(peaks[0][j],  precursorMass + massShift/3, 2) < 10){
							System.out.println("Scan Number " + scanNum + " has charge 3 isotope");
							return true;
						}

					}
		}
		return false;
	}
	
	private int getPrevScan1(int ScanNum){
		int scan = ScanNum - 1;
		while(scan > 0){
			Scan MS1 = this.parser.rap(scan);
			if(MS1.getHeader().getMsLevel() == 1){
				return scan;
			}
			scan--;
		}
		return -1;
	}
	
	public int getPrevScan(int ScanNum, int MsLevel){
		int scan = ScanNum - 1;
		while(scan > 0){
			Scan MS = this.parser.rap(scan);
			//System.out.println("scan number: " + scan);
			if(MS != null && MS.getHeader().getMsLevel() == MsLevel){
				return scan;
			}
			scan--;
		}
		return -1;
	}
	
	
	public double getPIF(int scanNum){
		Scan spectrum = this.parser.rap(scanNum);
		Scan MS1 = this.parser.rap(getPrevScan1(spectrum.getHeader().getNum()));
		double ppm = 20;
		double windowWidth = 2.0;
		int numIsoPeaks = 4;
		double precursor = spectrum.getHeader().getPrecursorMz();
		int charge = spectrum.getHeader().getPrecursorCharge();
		return getPIF(scanNum, precursor, ppm, -1*windowWidth, windowWidth, numIsoPeaks, charge)[2];
	}
	
	public double[] getPIF(int scanNum, double precursor, double ppm, double windowLeft, double windowRight, int numIsoPeaks, int charge){
		Scan spectrum = this.parser.rap(scanNum);
		Scan MS1 = this.parser.rap(getPrevScan1(spectrum.getHeader().getNum()));
		//System.out.println("MS1: " + MS1.getHeader().getNum());
		double[][] massIntList = MS1.getMassIntensityList();
		double totalInt = 0.000001;
		double precursorInt = 0.0;
		double[] precursorProfile = new double[numIsoPeaks+1];
		for(int i = -1; i < numIsoPeaks; i++){
			precursorProfile[i+1] = precursor+(1.0*i)/charge;
		}
		
		double max = 0.0;
		for(int i = 0; i < massIntList[1].length; i++){
			max = max > massIntList[1][i] ? max : massIntList[1][i];
		}
				
		for(int i = 0; i < massIntList[0].length; i++){
			double diff = massIntList[0][i] - precursor;
			//System.out.println("diff " + diff);
			if(diff > windowLeft && diff < windowRight){
				if(massIntList[1][i]/max > 0.005){
					totalInt += massIntList[1][i];
					for(int j = 0; j < precursorProfile.length; j++){
						double diff2 = massIntList[0][i] - precursorProfile[j];
						if(Math.abs(diff2)*1000000/precursor < ppm){
							precursorInt += massIntList[1][i];
						}
					}
				}
			}
		}
		//System.out.println("Scan:\t" + scanNum + "\tprecursor:\t" +  precursor + "\t" + charge 
		//		+ "\tInt ratio:\t" + precursorInt + "\t" + totalInt + "\t" + (precursorInt/totalInt));
		return new double[]{precursorInt, totalInt, precursorInt/totalInt};
	}
	
	@Override
	public boolean hasNext() {
		return !(this.current == null);	
	}

	@Override
	public Spectrum next() {
		Spectrum prev = this.current;
		currentScan += 1;
		for(; currentScan < this.maxScan; currentScan++){
			Scan scan = parser.rap(currentScan);
			//wait a little bit before reading another scan, in linux this sometime cause too many file handles to be open
			//System.out.println("reading spectrum at: " + currentScan);
			try{
				parser.wait(100);
			}catch(Exception e){
				
			}
			if(scan != null && scan.getHeader().getMsLevel() >= 2){
				this.current = getSpectrum(currentScan);
				//this.current.score = scan.getHeader().getPrecursorIntensity();
				double totalInt = this.current.sumMagnitude();
				//current.getHeader().getTotIonCurrent()
				if(scan.getHeader().getMsLevel() == 3){
					//Scan ms2 = parser.rap(currentScan-1);
					//this.current.charge = ms2.getHeader().getPrecursorCharge()-1;
				}
//			System.out.println("readed: " + scan.getHeader().getNum() + "\t"
//					    +"MS"+scan.getHeader().getMsLevel() + "\t"
//						+ scan.getHeader().getPrecursorMz() + "\t" 
//						+ scan.getHeader().getPrecursorCharge() + "\t"
//						+ current.charge 
//						+ "\t" + scan.getHeader().getPrecursorIntensity() 
//						+"\t" + totalInt +"\t" + scan.getHeader().getTotIonCurrent());
				return prev;
			}			
		}
		this.current = null;
		return prev;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Spectrum readSpectrumByIndex(int index) {
		return this.getSpectrum(index);
	}
	                                                
		

}
