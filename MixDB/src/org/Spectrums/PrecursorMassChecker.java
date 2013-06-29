package org.Spectrums;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import Utils.FileIOUtils;
import org.systemsbiology.jrap.*;
import org.systemsbiology.jrap.stax.*;


/**
 * Check precursor mass for MS/MS spectrum annotation
 * in the case of high-accuracy MS scan we can use this
 * as a first-level of validation of IDs returned by our search algorithm
 * Alternatively we can also use  this as a preprocessing step to select select
 * candidate for scoring, since precurosr recorded in MS/MS maybe not be 
 * corrected or the first isotopic mass of peptide being fragment we
 * wil scan the MS1 profiel for matches
 * @author Jian Wang
 *
 */
public class PrecursorMassChecker {
	private static Peak dummy = new Peak(0,0);
	private String filename;
	private MSXMLParser parser;
	private MZXMLReader reader;
	private Map<String, String> annotation;
	private List<String> annotations;
	private boolean printDetail = false;
	private SortedMap<Integer, SortedMap<Double, Peak>> spectraMap;
	private SortedMap<Integer, Spectrum> spectraMap2;
	private Map<Integer, Double> scanRTMapper;
	private double systematicError=0;
	
	public double getSystematicError() {
		return systematicError;
	}

	public void setSystematicError(double systematicError) {
		this.systematicError = systematicError;
	}

	public PrecursorMassChecker(String filename, String annotationFile){
		this.filename = filename;
		this.parser = new MSXMLParser(this.filename);
		loadAnnotation(annotationFile);
	}
	
	public PrecursorMassChecker(String filename){
		this.filename = filename;
		this.parser = new MSXMLParser(this.filename);
	}
	
	private void loadAnnotation(String annotation){
		//this.annotation = Utils.FileIOUtils.createTableFromFile(annotation, 0, 1);
		this.annotations = Utils.FileIOUtils.createListFromFile(annotation);
		System.out.println("read in " + this.annotations.size() + " annotation"); 
	}
	
	private String[] parseResultLine(String line){
		String[] tokens = line.split("\\s+");
		return new String[]{tokens[3], tokens[7]};
	}
	
	private String[] parseMixtureResultLine(String line){
		String[] tokens = line.split("\\s+");
		return new String[]{tokens[4], tokens[5], tokens[7]};
	}
	
	//convert all MS1 spectra into sortedmap so we can extract peaks 
	//in a more efficient way
	public SortedMap<Integer, SortedMap<Double, Peak>> generateMS1Map(){
		this.spectraMap = new TreeMap<Integer, SortedMap<Double, Peak>>();
		for(int i = 1; i < this.parser.getScanCount(); i++){
			Scan current = this.parser.rap(i);
			List<Peak> peakList;
			if(current.getHeader().getMsLevel() == 1){
				peakList = getPeakList(current);
				SpectrumUtil.toRelativeIntensity(peakList);
				removeNoisePeaks(peakList, 0.02);
				SortedMap<Double, Peak> spectraMap = new TreeMap<Double, Peak>();
				for(int j = 0; j < peakList.size(); j++){
					Peak p = peakList.get(j);
					spectraMap.put(p.getMass(), p);
				}
				this.spectraMap.put(i, spectraMap);
			}
		}
		return this.spectraMap;
	}
	
	public SortedMap<Integer, Spectrum> generateMS2Map(){
		this.spectraMap2 = new TreeMap<Integer, Spectrum>();
		for(int i = 1; i < this.parser.getScanCount(); i++){
			Scan current = this.parser.rap(i);
			List<Peak> peakList;
			if(current.getHeader().getMsLevel() == 2){
				Spectrum s = MZXMLReader.getSpectrum(current);
				this.spectraMap2.put(i, s);
			}
		}
		return this.spectraMap2;
	}
	
	public List<Peak> getPeakList(Scan s){
		List<Peak> pList = new ArrayList<Peak>();
		double[][] peaks = s.getMassIntensityList();
		for(int i = 0; i < peaks[0].length; i++){
			Peak p = new Peak(peaks[0][i], peaks[1][i]);
			pList.add(p);
		}
		//System.out.println("last peak from " + s.getNum() + "\t" + peaks[0][peaks[0].length-1] + "\t" + peaks[1][peaks[1].length-1]);
		return pList;
	}
	//for a particular MS2, check whether annotation give match the
	//precursors in the corresponding MS1 
	public int matchPrecursorProfile(int scan, Peptide p, double tolerance){
		if(this.spectraMap == null){
			this.spectraMap = this.generateMS1Map();
		}
		if(this.spectraMap2 == null){
			this.spectraMap2 = this.generateMS2Map();
		}
	//	Scan MS2 = this.parser.rap(scan);
		Spectrum MS2 = this.spectraMap2.get(scan);
		int scanMS1 = getPrecursorScan(scan);
		//Scan MS1 = parser.rap(scanMS1);
		SortedMap MS1 = this.spectraMap.get(scanMS1);
		//System.out.println("MS1 is null: " + scanMS1);
		double[] diff = computeMinDiff(p, MS1, 0, false, MS2.parentMass, 3);
		//double[] signedDiff = computeSignedMinDiff(MS1, 0, false, MS2.parentMass, 3);
		//System.out.println("Scan " + scan + ":\t" + p + "\t" + p.getParentmass() + "\tms1 error:\t" + diff[0] + "\t" + diff[1] + "\t" + diff[2] + "\t" + diff[3]);
		for(int i = 0; i < diff.length; i++){
			if(diff[i] > tolerance){
				//System.out.println("scan: " + scan + "has match: " + i);
				return i;
			}
		}
		return diff.length;
	}
	
	//given a peptide, check whether the particular precursor
	//show up in the run, and where does it show up
	public int matchPeptidePrecursorProfile(Peptide p, double tolerance){
		for(int i = 700; i < 2300; i++){
			Scan current = this.parser.rap(i);
			if(current.getHeader().getMsLevel() == 1){
				double[] diff = computeMinDiff(p, current, 0, false, p.getParentmass(), 3);
				int j = 0;
				for(j = 0; j < diff.length; j++){
					if(diff[j] > tolerance){
						//System.out.println("scan: " + scan + "has match: " + i);
						break;
					}
				}
				if(j >= 3){
					System.out.println("Peptide: " + p + "\t" + p.getCharge() + "\t" 
							+ p.getParentmass() + " show up in Scan:\t" + i + "\twith match peak\t" +  j);
				}
			}
		}
		return 0;
	}
	
	public int matchPeptidePrecursorProfile2(Peptide p, double tolerance){
		if(this.spectraMap == null){	
			this.generateMS1Map();
		}
		int matchCount=0;
		for(Iterator<Integer> it = this.spectraMap.keySet().iterator(); it.hasNext();){
			Integer scanNum = it.next();
			SortedMap spectrum = this.spectraMap.get(scanNum);		
			double[] diff = computeMinDiff(p, spectrum, 0, false, p.getParentmass(), 3);
			int j = 0;
			for(j = 0; j < diff.length; j++){
				if(diff[j] > tolerance){
					//System.out.println("scan: " + scan + "has match: " + i);
					break;
				}
			}
			if(j >= 3){
				System.out.println("Peptide: " + p + "\t" + p.getCharge() + "\t" 
						+ p.getParentmass() + " show up in Scan:\t" + scanNum + "\twith match peak\t" +  j);
				matchCount++;
			}
		}
		return matchCount;
	}
	
	public int[] matchPeptidePrecursorProfilePair(Peptide p, double tolerance, double isotopeOffset){
		List[] matchStat = matchPeptidePrecursorProfilePairDetail(p, tolerance, isotopeOffset);
		return new int[]{matchStat[0].size(), matchStat[1].size(), matchStat[2].size(),
						matchStat[3].size(), matchStat[4].size(), matchStat[5].size(), matchStat[6].size()};
	}
	
	public List[] matchPeptidePrecursorProfilePairDetail(Peptide p, double tolerance, double isotopeOffset){
		List [] matchStat = new List[7];
		//initialize return list;
		for(int i = 0; i < matchStat.length; i++){
			matchStat[i] = new ArrayList();
		}
		
		if(this.spectraMap == null){	
			this.generateMS1Map();
		}
		
		if(this.spectraMap2 == null){
			this.generateMS2Map();
		}
		
		if(this.scanRTMapper == null){
			if(this.reader == null){
				this.reader = new MZXMLReader(this.parser);
			}
			this.scanRTMapper = reader.getRTScanMapping();
		}
		
		for(Iterator<Integer> it = this.spectraMap.keySet().iterator(); it.hasNext();){
			Integer scanNum = it.next();
			SortedMap spectrum = this.spectraMap.get(scanNum);
			Object[] ret1 = computeMinDiff(spectrum, 0, false, p.getParentmass(), 3, p.getCharge());
			Object[] ret2 = computeMinDiff(spectrum, 0, false, p.getParentmass()+(isotopeOffset/p.getCharge()), 3, p.getCharge());
			//Object[] ret2 = ret1;
			double[] diff = (double[])ret1[0];
			double[] diff2 = (double[])ret2[0];
			Peak[] peaks1 = (Peak[])ret1[1];
			Peak[] peaks2 = (Peak[])ret2[1];
//			System.out.println("Scan: " + scanNum + " diff: " + diff[0] + "\t" + diff[1] + "\t" + diff[2]);
			int j = 0;
			for(j = 0; j < diff.length; j++){
				if(diff[j] > tolerance){
					//System.out.println("scan: " + scan + "has match: " + i);
					break;
				}
			}
			int k = 0;
			for(k = 0; k < diff2.length; k++){
				if(diff2[k] > tolerance){
					//System.out.println("scan: " + scan + "has match: " + i);
					break;
				}
			}
			
			int matchType1 = 0, matchType2 = 0;
			if(j > 2){
				System.out.println("checking scan: " + scanNum.intValue());
				matchType1 = checkMS2(scanNum.intValue(), p);
				if(matchType1 >= 1){
					matchStat[3].add(scanNum);
				}
				if(matchType1 == 2){
					matchStat[4].add(scanNum);
				}
				matchStat[0].add(scanNum);
			}
			
			if(k>2){
				Peptide isoD12 = new Peptide(p);
				if(isoD12.getPtmmasses().length > 0){
					isoD12.getPtmmasses()[0] += (isotopeOffset);
					isoD12.setParentmass(isoD12.getParentmass() + (isotopeOffset/p.getCharge()));
					matchType2 = checkMS2(scanNum.intValue(), isoD12);
				}
				if(p instanceof LinkedPeptide){
					LinkedPeptide lp2 = new LinkedPeptide((LinkedPeptide)p);
					Peptide p1_D12 = lp2.peptides[0];
					Peptide p2_D12 = lp2.peptides[1];
					p1_D12.getPtmmasses()[0] += (isotopeOffset);
					p1_D12.setParentmass(isoD12.getParentmass() + isotopeOffset);
					p2_D12.getPtmmasses()[0] += (isotopeOffset);
					p2_D12.setParentmass(isoD12.getParentmass() + isotopeOffset);
					lp2.setParentmass(isoD12.getParentmass() + (isotopeOffset/p.getCharge()));
					matchType2 = checkMS2(scanNum.intValue(), lp2);
				}	
				
				if(matchType2 >= 1){
					matchStat[5].add(scanNum);
				}
				if(matchType2 == 2){
					matchStat[6].add(scanNum);
				}
				matchStat[1].add(scanNum);
			}
			if(j > 2 && k > 2){ //both isotopic pairs show up
				matchStat[2].add(scanNum);
			}
			if(j > 2 || k > 2){
				diff = computeSignedMinDiff(spectrum, 0, false, p.getParentmass(), 3, p.getCharge());
				diff2 = computeSignedMinDiff(spectrum, 0, false, p.getParentmass()+(isotopeOffset/p.getCharge()), 3, p.getCharge());
				double rttime = this.scanRTMapper.get(scanNum);
				System.out.println("Peptide: " + p + "\t" + p.getCharge() + "\t" 
						+ p.getParentmass() + " show up in Scan:\t" + scanNum + "\ttime:\t" + rttime + "\twith match peak\t" +  j 
						+ "\twith heavy isotopic peaks:\t" + k  + "\tMS2 stat: " + matchType1 +"\t" + matchType2
						+ "\tMS1 error: " + diff[0] + "\t" + diff[1] + "\t" + diff[2] + "\t" + diff[3] + "\t" + diff[4] 
						+ "\tMS1 error: " + diff2[0] + "\t" + diff2[1] + "\t" + diff2[2] + "\t" + diff2[3] + "\t" + diff2[4]
						+ "\tMS1 Intensity: " + peaks1[0].getIntensity() + "\t" + peaks1[1].getIntensity() + "\t" + peaks1[2].getIntensity()
						+ "\tMS1 Intensity: " + peaks2[0].getIntensity() + "\t" + peaks2[1].getIntensity() + "\t" + peaks2[2].getIntensity() +"\t"
						+ "\tMS1 Mass: " + peaks1[0].getMass() + "\t" + peaks1[1].getMass() + "\t" + peaks1[2].getMass() + "\t"
						+ "\tMS1 Mass: " + peaks2[0].getMass() + "\t" + peaks2[1].getMass() + "\t" + peaks2[2].getMass() +"\t"
						);
			}
			
		}
		System.out.println("both show up in scan counts: " + matchStat[2].size());
		return matchStat;
	}
	
	//check MS2 match-info, return match type
	public int checkMS2(int MS1Scan, Peptide p){
		int matchType = 0;
		List<Spectrum> MS2List = getAssociatedMS2(MS1Scan);
		for(int i = 0; i < MS2List.size(); i++){
			matchType = 0;
			Spectrum s = MS2List.get(i);
			s.windowFilterPeaks(10, 25);
			s.computePeakRank();
			s.charge = p.getCharge();
			if(Math.abs(s.parentMass - p.getParentmass()) < 1){ 
					//Math.abs(s.parentMass - (Mass.C13-Mass.C12)/s.charge - p.getParentmass()) < 0.05){
				matchType++;
			}else{
				continue;
			}
			if(p instanceof LinkedPeptide){
				Peptide p1 = ((LinkedPeptide)p).peptides[0];
				Peptide p2 = ((LinkedPeptide)p).peptides[1];
				String[] peps = p.getPeptide().split("--");
				System.out.println("peptide1 is: " + peps[0] + "\t" + "peptides2 is: " + peps[1]);
				TheoreticalSpectrum linkedSpect = new TheoreticalSpectrum(p1, p2, p.getCharge(), false); //note linkedSpec parent mass seems not correct in this case
				double[] stat = linkedSpect.analyzeMixtureAnnotation(s, peps[0], peps[1], 0.05);
				System.out.println("Spectrum: " + s.spectrumName+ "\t"  + s.parentMass + " has best match: " +  p1 + "--" + p2  + "\t" +  p.getParentmass()
						+ "\t" + stat[0] + "\t" + stat[1] + "\t"
						+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
						+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
						+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12]);
				if(stat[0] > 0.5 && (stat[1] + stat[2] > 0.85 && stat[3] + stat[4] > 0.85 && matchType > 0)){
					matchType++;
				}
			}else{
				TheoreticalSpectrum th = new TheoreticalSpectrum(p);
				double[] stat = th.analyzeAnnotation(s, p.getPeptide(), 0.3);
				System.out.println("Spectrum: " + s.spectrumName+ "\t"  + s.parentMass + " has best match: " +  p + "\t" +  p.getParentmass()
						+ "\t" + stat[0] + "\t" + stat[1] + "\t"
						+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
						+ stat[5]);
				if(stat[0] > 0.6 && stat[1] + stat[2] > 1 && matchType > 0){
					matchType++;
				}
			}
			if(matchType > 1){
				return matchType;
			}
		}
		return matchType;
	}
	
	//get all the MS2 associate with this MS1
	public List<Spectrum> getAssociatedMS2(int MS1Scan){
		List<Spectrum> specList = new ArrayList<Spectrum>();
		for(int i = MS1Scan+1; i < this.parser.getScanCount(); i++){
			Scan current = this.parser.rap(i);
			if(current.getHeader().getMsLevel() == 2){
				if(this.spectraMap2.containsKey(current.getHeader().getNum())){
					specList.add(this.spectraMap2.get(current.getHeader().getNum()));
				}
			}else{
				break;
			}
		}
		return specList;
	}
	
	public void checkPrecursorMZ(){
		for(Iterator<String> it = annotations.iterator(); it.hasNext();){
			String line = it.next();
			String[] match = parseResultLine(line);
			int scanNum = Integer.parseInt(match[0]);//.split("_")[1]);
			Scan current = parser.rap(scanNum);
			Scan currentMS1 = parser.rap(getPrecursorScan(scanNum));
			int charge=0;
			//System.out.println("match is : " + match[1]);
			Peptide p;
			if(match[1].split("\\.").length == 1){ //no charge info
				charge = getCharge(current, match[1]);
				System.out.println("guessed charge: " + charge);
				p = new Peptide(match[1]+"."+charge);
			}else{
				p = new Peptide(match[1]);
			}
			if(this.printDetail){
				System.out.print(line);
				System.out.print("\t");
			}
			computeMinDiff(p, scanNum, 0, true);
//			if(p.getPeptide().contains("N") || p.getPeptide().contains("Q")){
//				if(this.printDetail){
//					System.out.print(line);
//					System.out.print("\t");
//				}
//				p.insertPTM(1, 0.9840);
//				computeMinDiff(p, scanNum, 0, true);
//			}
		}
	}	
	
	
	public void checkLinkedPrecursorMZ(){
		for(Iterator<String> it = annotations.iterator(); it.hasNext();){
			String line = it.next();
			String[] match = parseResultLine(line);
			int scanNum = Integer.parseInt(match[0]);//.split("_")[1]);
			Scan current = parser.rap(scanNum);
			Scan currentMS1 = parser.rap(getPrecursorScan(scanNum));
			int charge=0;
			//System.out.println("match is : " + match[1]);
			Peptide p;
			match[1] = match[1].replaceAll("[0-9\\.\\+]", "");
			if(match[1].split("\\.").length == 1){ //no charge info
				charge = getCharge(current, match[1]);
				System.out.println("guessed charge: " + charge);
				p = new LinkedPeptide(match[1], current.getHeader().getPrecursorCharge());
			}else{
				p = new Peptide(match[1], current.getHeader().getPrecursorCharge());
			}
			if(this.printDetail){
				System.out.print(line);
				System.out.print("\t");
			}
			computeMinDiff(p, scanNum, 0, true);
//			if(p.getPeptide().contains("N") || p.getPeptide().contains("Q")){
//				if(this.printDetail){
//					System.out.print(line);
//					System.out.print("\t");
//				}
//				p.insertPTM(1, 0.9840);
//				computeMinDiff(p, scanNum, 0, true);
//			}
		}
	}	
	
	public int getCharge(Scan scan, String p){
		int min = 1, max = 6;
		Peptide[] peps = new Peptide[max-min+1];
		for(int i = min; i < max+1; i++){
			peps[i-min] = new Peptide(p, i);
		}
		double minDff=1000000;
		int minInd=0;
		for(int i = min; i < max+1; i++){
			double diff = Math.abs(scan.getHeader().getPrecursorMz()-peps[i-min].getParentmass());
			//System.out.println("diff is: " + diff);
			minInd = diff < minDff ? i:minInd;
			minDff = diff < minDff ? diff:minDff;
		}
		System.out.println("mass difference is: " + minDff);
		return minInd;
	}
	
	public void checkMixturePrecursorMZ(){
		for(Iterator<String> it = annotations.iterator(); it.hasNext();){
			String line = it.next();
			String[] match = parseMixtureResultLine(line);
			int scanNum = Integer.parseInt(match[0]);
			Scan current = parser.rap(scanNum);
			Scan currentMS1 = parser.rap(getPrecursorScan(scanNum));
			String[] peptides = {match[1], match[2]};
			
			for(int i = 0; i < peptides.length; i++){
				//System.out.println("peptides is: " + peptides[i]);
				if(this.printDetail){
					System.out.print(line);
					System.out.print("\t");	
				}
				Peptide p = new Peptide(peptides[i]);
				computeMinDiff(p, scanNum, 0, true);
				//System.out.println("peptide become: " + p.getPeptide());
			}
		}
	}
		
	public boolean isPrintDetail() {
		return printDetail;
	}


	public void setPrintDetail(boolean printDetail) {
		this.printDetail = printDetail;
	}

	
	private double[] computeMinDiff(Peptide p, int scanNum, int mode){
		return computeMinDiff(p, scanNum, mode, false);
	}
	
	private double[] computeMinDiff(Peptide p, int scanNum, int mode, boolean detail){
		double[] profile = getProfile(p);
		int scanMS1 = getPrecursorScan(scanNum);
		Scan MS1 = parser.rap(scanMS1);
		Scan MS2 = parser.rap(scanNum);
		double[] closest = new double[profile.length];
		double[] minDiff = new double[profile.length];
		for(int i = 0; i < profile.length; i++){
			closest[i] = 10000000;
			minDiff[i] = 10000000;
		}
		double[][] peaks = MS1.getMassIntensityList();
		for(int i = 0; i < peaks[0].length; i++){
			for(int j = 0; j < profile.length; j++){
				//double diff = Math.abs(diff(peaks[0][i],  profile[j], 2));
				double diff = signedDiff(peaks[0][i],  profile[j], 2);
				closest[j] = Math.abs(diff) < Math.abs(minDiff[j]) ? peaks[0][i] : closest[j];
				minDiff[j] = Math.abs(diff) < Math.abs(minDiff[j]) ? diff(peaks[0][i],  profile[j], 2) : minDiff[j];
			}
		}
		
		if(detail){
			System.out.print("MS2 scan " + scanNum + "\t" +  MS2.getHeader().getPrecursorMz() + "\t" + MS2.getHeader().getPrecursorCharge() + " MS1 scan " + getPrecursorScan(scanNum) 
				+ " annotation: " + p.getPeptide() + "\t" 
				+	"closest match ");
			int count = 0;
			System.out.print("profile:\t");
			for(int i = 0; i < profile.length; i++){
				System.out.print(profile[i] + "\t");
			}
			System.out.print("closest:\t");
			for(int i = 0; i < closest.length; i++){
				System.out.print(closest[i] + "\t");
			}
			System.out.print("diff: " );
			for(int i = 0; i < minDiff.length; i++){
				if(Math.abs(minDiff[i]) < 5 && count==i){
					count++;
				}
				System.out.print(minDiff[i] + "\t");
			}
			System.out.print(" total matched: " + count + "\t");
		}
		int prevScan = getPrevScan1(MS1.getHeader().getNum());
		computeMinDiff(p, this.parser.rap(prevScan), 1, detail);
		int nextScan = getNextScan1(MS1.getHeader().getNum());	
		computeMinDiff(p, this.parser.rap(nextScan), 1, detail);
		if(detail){
			System.out.println();
		}
		return minDiff;

	}
	
	public void checkPrecurosrMass(Peptide p, int scanNum, int mode){
		double[] diff = computeMinDiff(p, scanNum, mode);
		int scanMS1 = getPrecursorScan(scanNum);
		Scan MS1 = parser.rap(scanMS1);
		
		
	}
	
	
	private double diff(double mass1, double mass2, int mode){
		if(mode == 1){
			return Math.abs(mass1 - mass2);
		}else {
			return Math.abs(((mass1 - mass2)*1000000 / mass2) - this.systematicError); 
		}

	}
	
	private double signedDiff(double mass1, double mass2, int mode){
		if(mode == 1){
			return mass1 - mass2;
		}else {
			return ((mass1 - mass2)*1000000 / mass2) - this.systematicError; 
		}

	}
	
	private double[] computeMinDiff(Peptide p, Scan MS1, int mode, boolean detail){
		double[] profile = getProfile(p);
		double[] closest = new double[profile.length];
		double[] minDiff = new double[profile.length];
		for(int i = 0; i < profile.length; i++){
			closest[i] = 10000000;
			minDiff[i] = 10000000;
		}
		double[][] peaks = MS1.getMassIntensityList();
		for(int i = 0; i < peaks[0].length; i++){
			for(int j = 0; j < profile.length; j++){
				//double diff = Math.abs(diff(peaks[0][i],  profile[j], 2));
				double diff = signedDiff(peaks[0][i],  profile[j], 2);
				closest[j] = Math.abs(diff) < Math.abs(minDiff[j]) ? peaks[0][i] : closest[j];
				minDiff[j] = Math.abs(diff) < Math.abs(minDiff[j]) ? diff(peaks[0][i],  profile[j], 2) : minDiff[j];
			}
		}
//		System.out.print("MS scan " + scanNum + " MS1 scan " + getPrecursorScan(scanNum) 
//				+ " annotation: " + p.getPeptide() + "\t" 
//				+	"closest match ");
		if(detail){
			System.out.print("MS1 scan: " + MS1.getHeader().getNum()+ "\t");
			int count = 0;
			System.out.print("profile:\t");
			for(int i = 0; i < profile.length; i++){
				System.out.print(profile[i] + "\t");
			}
			System.out.print("closest:\t");
			for(int i = 0; i < closest.length; i++){
				System.out.print(closest[i] + "\t");
			}
			System.out.print("diff:\t");
			for(int i = 0; i < minDiff.length; i++){
				if(Math.abs(minDiff[i]) < 5 && count==i){  
					count++;
				}
				System.out.print(minDiff[i] + "\t");
			}
			System.out.print(" total matched: " + count + "\t");
		}
		return minDiff;
	}

	private double[] computeMinDiff(Peptide p, Scan MS1, int mode, boolean detail, double center, double width){
		double[] profile = getProfile(p);
		double[] closest = new double[profile.length];
		double[] minDiff = new double[profile.length];
		for(int i = 0; i < profile.length; i++){
			closest[i] = 10000000;
			minDiff[i] = 10000000;
		}	
		List<Peak> peakList = extractLocalPeaks(MS1, 0, 5000); //extract all peaks
		removeNoisePeaks(peakList, 0.02);
		List<Peak> localPeakList = extractLocalPeaks(peakList, center, width);
		removeNoisePeaks(localPeakList, 0.05);			
		for(int i = 0; i < localPeakList.size(); i++){
			for(int j = 0; j < profile.length; j++){
				double currentmass =localPeakList.get(i).getMass();
				double diff = Math.abs(diff(currentmass,  profile[j], 2));
				closest[j] = diff < minDiff[j] ?  currentmass : closest[j];
				minDiff[j] = diff < minDiff[j] ? diff : minDiff[j];
			}
		}
		return minDiff;
	}
	
	private double[] computeSignedMinDiff(Peptide p, Scan MS1, int mode, boolean detail, double center, double width){
		double[] profile = getProfile(p);
		double[] closest = new double[profile.length];
		double[] minDiff = new double[profile.length];
		for(int i = 0; i < profile.length; i++){
			closest[i] = 10000000;
			minDiff[i] = 10000000;
		}	
		List<Peak> peakList = extractLocalPeaks(MS1, 0, 5000); //extract all peaks
		removeNoisePeaks(peakList, 0.02);
		List<Peak> localPeakList = extractLocalPeaks(peakList, center, width);
		removeNoisePeaks(localPeakList, 0.05);			
		for(int i = 0; i < localPeakList.size(); i++){
			for(int j = 0; j < profile.length; j++){
				double currentmass =localPeakList.get(i).getMass();
				double diff = signedDiff(currentmass,  profile[j], 2);
				closest[j] = Math.abs(diff) < Math.abs(minDiff[j]) ?  currentmass : closest[j];
				minDiff[j] = Math.abs(diff) < Math.abs(minDiff[j]) ? diff : minDiff[j];
			}
		}
		return minDiff;
	}
	
	private double[] computeMinDiff(Peptide p, SortedMap scan, int mode, boolean detail, double center, double width){
		double[] profile = getProfile(p);
		Peak[] closest = new Peak[profile.length];
		double[] minDiff = new double[profile.length];
		for(int i = 0; i < profile.length; i++){
			minDiff[i] = 10000000;
		}
		List<Peak> localPeakList = new ArrayList<Peak>(); //extract all peaks
		localPeakList.addAll(scan.subMap(center-width, center+width).values());
//		System.out.println("interaval: " + (center-width) + " :" + (center + width));
//		for(int i = 0; i < localPeakList.size(); i++){
//			System.out.print(localPeakList.get(i) + "\t");
//		}
//		System.out.println();
		removeNoisePeaks(localPeakList, 0.05);			
		for(int i = 0; i < localPeakList.size(); i++){
			for(int j = 0; j < profile.length; j++){
				Peak current = localPeakList.get(i);
				double currentmass =localPeakList.get(i).getMass();
				double diff = Math.abs(diff(currentmass,  profile[j], 2));
				closest[j] = diff < minDiff[j] ?  current : closest[j];
				minDiff[j] = diff < minDiff[j] ? diff : minDiff[j];
			}
		}
//		if(minDiff[0] != 10000000 && closest[0].getIntensity() < 0.15){
//			minDiff[0] = 1000;
//		}
//		for(int i = 0; i < closest.length; i++){
//			System.out.print(closest[i] + "\t");
//		}
//		System.out.println();
		return minDiff;
	}
	
	private Object[] computeMinDiff(SortedMap scan, int mode, boolean detail, double center, double width, int charge){
		double[] profile = getProfile(center, charge);
		Peak[] closest = new Peak[profile.length];
		double[] minDiff = new double[profile.length];
		for(int i = 0; i < profile.length; i++){
			minDiff[i] = 10000000;
			closest[i] = dummy;
		}
		List<Peak> localPeakList = new ArrayList<Peak>(); //extract all peaks
		localPeakList.addAll(scan.subMap(center-width, center+width).values());
		removeNoisePeaks(localPeakList, 0.05);
		for(int i = 0; i < localPeakList.size(); i++){
			for(int j = 0; j < profile.length; j++){
				Peak current = localPeakList.get(i);
				double currentmass =localPeakList.get(i).getMass();
				double diff = Math.abs(diff(currentmass,  profile[j], 2));
				closest[j] = diff < minDiff[j] ?  current : closest[j];
				minDiff[j] = diff < minDiff[j] ? diff : minDiff[j];
			}
		}
//		if(minDiff[0] != 10000000 && closest[0].getIntensity() < 0.10){
//			minDiff[0] = 1000;
//		}
		
		if(minDiff[0] != 10000000){
			if(checkProfileIntensity(closest)){
				return new Object[]{minDiff, closest};
			}else{
				minDiff[0] = 1000;
			}
		}
		return new Object[]{minDiff, closest};
	}
	
	private double[] computeSignedMinDiff(SortedMap scan, int mode, boolean detail, double center, double width, int charge){
		double[] profile = getProfile(center, charge);
		Peak[] closest = new Peak[profile.length];
		double[] minDiff = new double[profile.length];
		for(int i = 0; i < profile.length; i++){
			minDiff[i] = 10000000;
		}
		List<Peak> localPeakList = new ArrayList<Peak>(); //extract all peaks
		localPeakList.addAll(scan.subMap(center-width, center+width).values());
		removeNoisePeaks(localPeakList, 0.05);			
		for(int i = 0; i < localPeakList.size(); i++){
			for(int j = 0; j < profile.length; j++){
				Peak current = localPeakList.get(i);
				double currentmass =localPeakList.get(i).getMass();
				double diff = signedDiff(currentmass,  profile[j], 2);
				closest[j] = Math.abs(diff) < Math.abs(minDiff[j]) ?  current : closest[j];
				minDiff[j] = Math.abs(diff) < Math.abs(minDiff[j]) ? diff : minDiff[j];
			}
		}
//		if(minDiff[0] != 10000000 && closest[0].getIntensity() < 0.10){
//			minDiff[0] = 1000;
//		}
		
		if(minDiff[0] != 10000000){
			if(checkProfileIntensity(closest)){
				return minDiff;
			}else{
				minDiff[0] = 1000;
			}
		}
		return minDiff;
	}
	
	private boolean checkProfileIntensity(Peak[] closest){
		boolean test1 = closest[0].getIntensity() > closest[2].getIntensity();
		boolean test2 = closest[1].getIntensity() > closest[2].getIntensity();
		boolean test3 = closest[0].getIntensity() > closest[1].getIntensity();
		boolean test4 = closest[1].getIntensity() / closest[0].getIntensity() < 1.5;
		//return test1 && test2 && (test3 || test4);
		//return test3;
		return true;
	}
	
	//extract peaks width Da away from the center mass
	private List<Peak> extractLocalPeaks(Scan MS1, double center, double width){
		List<Peak> peakList = new ArrayList<Peak>();
		double[][] peaks = MS1.getMassIntensityList();
		for(int i = 0; i < peaks[0].length; i++){
			if(Math.abs(peaks[0][i] - center) < width){
				Peak p = new Peak((double)peaks[0][i], (double)peaks[1][i]);
				peakList.add(p);
			}
		}	
		return peakList;
	}
	
	private List<Peak> extractLocalPeaks(List<Peak> peakList, double center, double width){
		List<Peak> localPeakList = new ArrayList<Peak>();
		for(int i = 0; i < peakList.size(); i++){
			Peak currentPeak = peakList.get(i);
			if(Math.abs(currentPeak.getMass() - center) < width){
				localPeakList.add(currentPeak);
			}
		}	
		return peakList;
	}
	//remove peaks the is not at X% of the highest intensity peaks in the list
	private void removeNoisePeaks(List<Peak> peakList, double threshold){
		double maxIntensity=0;
		for(int i = 0; i < peakList.size(); i++){
			double intensity = peakList.get(i).getIntensity();
			maxIntensity = maxIntensity > intensity ? maxIntensity : intensity;
		}
		//System.out.println("max is: " + maxIntensity);
		for(Iterator<Peak> it = peakList.iterator(); it.hasNext();){
			Peak current = it.next();
			if(current.getIntensity() / maxIntensity < threshold){
				//System.out.println("fraction is: " + current.getIntensity() / maxIntensity);
				//	System.out.println("removining " + current);
				it.remove();
			}
		}
	}
	
//	private void removeNoisePeaks(List<Peak> peakList, double threshold){
//		//System.out.println("max is: " + maxIntensity);
//		for(Iterator<Peak> it = peakList.iterator(); it.hasNext();){
//			Peak current = it.next();
//			if(current.getIntensity()  < threshold){
//				it.remove();
//			}
//		}
//	}
	
	private double[] getProfile(Peptide p){
		return getProfile(p.getParentmass(), p.getCharge());
	}

	private double[] getProfile(double mass, int charge){
		double[] profile = new double[5];
		for(int i = 0; i < profile.length; i++){
			profile[i] = mass + (i*1.00335/(double)charge);
		}
		return profile;
	}

	//this get the largest MS1 scan that preceeds the corresponding MS2 scan
	//for lack of better knowledge we assume scan start at one
	private int getPrecursorScan(int ScanNum){
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
	
	private int getNextScan1(int ScanNum){
		int scan = ScanNum + 1;
		while(scan < this.parser.getScanCount()+1){
			Scan MS1 = this.parser.rap(scan);
			if(MS1.getHeader().getMsLevel() == 1){
				return scan;
			}
			scan++;
		}
		return -1;
	}
	
	
	public static void computePrecursorMZError(String annotationFile){
		List<String> spectrumMass = 
			FileIOUtils.createListFromFile(annotationFile);
		//System.out.println("read in " + spectrumMass.size() + " lines");
		for(Iterator<String> idIter = spectrumMass.iterator(); idIter.hasNext();){
			String line = idIter.next();
			if(line.contains("#SpectrumFile")){
				line = idIter.next(); //skipping header
			}
			String[] tokens = line.split("\\s+");
			//String[] tokens = line.split("\t");
			//double actualMass = Double.parseDouble(tokens[20]);
			double actualMass = Double.parseDouble(tokens[8]);

			//System.out.println("actual mass is: " + actualMass);
			if(tokens[4].equals("0")){
				tokens[4] = "2";        //have no idea why inspect output charge 0, we replace with 2
			}
			//String pep = tokens[2].split("\\.")[1]+"."+tokens[4];
			String pep = tokens[5];
			//System.out.println("pep is: " + pep);
			
			Peptide p = new Peptide(pep);
			double[] profile = new double[]{p.getParentmass(), p.getParentmass()+(1/(double)p.getCharge()),
					p.getParentmass()+(2/(double)p.getCharge())};
			//double massError = actualMass - p.getParentmass();//computeMinDiff(actualMass, profile, 1);
			//double ppm = massError*1000000/actualMass;//computeMinDiff(actualMass, profile, 2);
			double massError = computeMinDiff(actualMass, profile, 1);
			double ppm = massError*1000000/actualMass;//computeMinDiff(actualMass, profile, 2);

			System.out.println(line + "\tprecursor-error:\t" + massError + "\t" + ppm);
		}
		
	}
		
	private static double computeMinDiff(double actual, double[] profile, int mode){
		double min = 100000;
		for(int i = 0; i < profile.length; i++){
			double diff = actual - profile[i];
			min = Math.abs(diff) < Math.abs(min) ? diff : min;
		}
		if(mode == 1){
			return min;
		}else{
			return min*1000000/actual;
		}
	}
	
	
	public static void testCheckPrecursorProfile(){
		String annotation = "..\\mixture_linked\\testAnnotation.txt";
		String fileName = "..\\mixture_linked\\human_heck_data/data/090121_NM_Trypsin_28.mzXML";
		PrecursorMassChecker checker = new PrecursorMassChecker(fileName, annotation);
		checker.setPrintDetail(true);
		checker.checkPrecursorMZ();
		//checker.checkMixturePrecursorMZ();
	}
	
	public static void testCheckLinkedPrecursorProfile(){
		String annotation = "..\\mixture_linked\\testAnnotation.txt";
		String fileName = "..\\mixture_linked\\linked_peptide_library\\sumo_lib\\20101008_Sumo_Library_4349_Bo.mzXML";
		PrecursorMassChecker checker = new PrecursorMassChecker(fileName, annotation);
		checker.setPrintDetail(true);
		checker.checkLinkedPrecursorMZ();
		//checker.checkMixturePrecursorMZ();
	}
	
	public static void main(String[] args){	
		testCheckPrecursorProfile();
		//testCheckLinkedPrecursorProfile();
		//computePrecursorMZError("..\\mixture_linked\\t");
	}
	
	


}
