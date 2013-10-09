package org.Spectrums;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import Utils.FileIOUtils;

public class HarklorParser {
	private String HarklorPath;
	public int minCharge = 2;
	public int maxCharge = 5;
	public int topByInt = 12;
	public int chargeInd = 2;
	public int mzInd = 4;
	
	public Map<Integer, SortedMap<Double, Double[]>> precursorList;  //we stored the precursors by  m/z for fast lookup
	public HarklorParser(String harklorFile){
		this.HarklorPath = harklorFile;
		this.precursorList = parseHarklorResults();
	}
	
	private Map<Integer, SortedMap<Double, Double[]>> parseHarklorResults(){
		Map<Integer, SortedMap<Double, Double[]>> harklorResults = new HashMap();
		List<String> results = FileIOUtils.createListFromFile(this.HarklorPath);
		int currentMS1 = 0;
		SortedMap<Double, Double[]> currPrecursorList =  null;
		for(int i = 0; i < results.size(); i++){
			String[] tokens = results.get(i).split("\\s+");
			//System.out.println("line is : " + results.get(i));
			if(tokens[0].equals("S")){
				if(currPrecursorList != null){
					//System.out.println("adding scan: " + currentMS1 + "\t" + currPrecursorList.size());
					harklorResults.put(currentMS1, currPrecursorList);
				}
				currentMS1 = Integer.parseInt(tokens[1]);
				currPrecursorList = new TreeMap<Double, Double[]>();
			}
			if(tokens[0].equals("P")){
				double mz = Double.parseDouble(tokens[mzInd]);
				double charge = Double.parseDouble(tokens[chargeInd]);
				double score = Double.parseDouble(tokens[8]);   //score
				double intensity = Double.parseDouble(tokens[3]); //base peak intensity
				if(charge >= this.minCharge && charge <= this.maxCharge){
					currPrecursorList.put(mz, new Double[]{mz, charge, score, intensity});
				}
			}
		}
		return harklorResults;
	}
	
	public static void getSWATHPrecursors(String hardklorPath, String swathFile){
		HarklorParser parser = new HarklorParser(hardklorPath);
		MZXMLReader reader = new MZXMLReader(swathFile);
		while(reader.hasNext()){
			Spectrum s = reader.next();
			int MS1Ind = SWATHUtils.getSWATHMS1Scan(s);
			if(parser.precursorList.containsKey(MS1Ind)){
				SortedMap<Double, Double[]> precursorList = parser.precursorList.get(MS1Ind);
				SortedMap<Double, Double[]> precursors = precursorList.subMap(s.parentMass - 5, s.parentMass + 26);
				for(Iterator<Double[]> it = precursors.values().iterator(); it.hasNext();){
					Double[] precursor = it.next();
					System.out.println(s.scanNumber + "\t" + precursor[0] + "\t" + precursor[1]);
				}
				//System.out.println("Swath-scans: " + s.scanNumber + "\thas detectable precursors: " + precursors.size());
			}
		}
	}
	
	//check if search results has harklor precursor present in the survey scans
	
	public static void checkHarklorPrecursor(String hardklorPath, String searchResults, String swathFile){
		checkHarklorPrecursor(hardklorPath, searchResults, swathFile, true);
	}
	
	public static void checkHarklorPrecursor(String hardklorPath, String searchResults, String swathFile, boolean checkCharge){
		HarklorParser parser = new HarklorParser(hardklorPath);
		List<String> results = Utils.FileIOUtils.createListFromFile(searchResults);
		MZXMLReader reader = new MZXMLReader(swathFile);
		double tolerance = 0.05;
		int matched = 0;
		for(int i = 0; i < results.size(); i++){
			String result = results.get(i);
			String[] tokens = result.split("\\t");
			if(tokens.length < 15 || result.startsWith("#")){
				continue;
			}
			int scanNumber = Integer.parseInt(tokens[1]);
			int MS1 = reader.getPrevScan(scanNumber, 1);
			System.out.println("MS1 " + MS1);
			String found = "NO-MATCH";
			if(parser.precursorList.containsKey(MS1)){
				double precursor = Double.parseDouble(tokens[5]);
				SortedMap<Double, Double[]> precursorList = parser.precursorList.get(MS1);
				//System.out.println("precursorlist size: " + precursorList.values().size());
				int charge = Integer.parseInt(tokens[6]);
				//System.out.println("MS1 is: " + MS1 + "\t" + precursorList.size());
				//int type = SWATHUtils.DALTON;
				//System.out.println("precursor: " + precursor);
				SortedMap<Double, Double[]> precursors = precursorList.subMap(precursor - tolerance, precursor + tolerance);
				if(precursors.size() > 0){
					for(Iterator<Double[]> it = precursors.values().iterator(); it.hasNext();){
						Double[] p = it.next();
						if(p[1] == charge || !checkCharge){
							found = "Match";
							break;
						}
					}
				}
				double precursorLeft = precursor-(Mass.C13-Mass.C12)/charge;
				//System.out.println("precursor: " + precursorLeft);
				precursors = precursorList.subMap(precursorLeft - tolerance, precursorLeft + tolerance);
				if(precursors.size() > 0){
					for(Iterator<Double[]> it = precursors.values().iterator(); it.hasNext();){
						Double[] p = it.next();
						if(p[1] == charge || !checkCharge){
							found = "Match";
							break;
						}
					}
				}
				double precursorRight = precursor+(Mass.C13-Mass.C12)/charge;
				//System.out.println("precursor: " + precursorRight);
				precursors = precursorList.subMap(precursorRight - 1, precursorRight + 1);
				//System.out.println("submap size: " + precursors.size());
				if(precursors.size() > 0){
					for(Iterator<Double[]> it = precursors.values().iterator(); it.hasNext();){
						Double[] p = it.next();
						if(p[1] == charge){
							found = "Match";
							break;
						}
					}
				}
				if(found.equals("MATCH")){
					matched++;
				}
				System.out.println(result+"\t"+found);
			}
		}
		System.out.println("found: " + matched +"/" + results.size() + ":\t" + (matched/(double)results.size()));
	}
	
	
	public static void main(String[] args){
		String hardklorPath = "../Hardklor/win32/14344_UPS1Ecolilysate_SWATH_5600_0.7min.hkout";
		String searchPath = "..//mixture_linked/ACG_swathdevelopment_14344_Swath_SSMs.txt";
		String swathPath = "..//mixture_linked/msdata/UPS_Ecoli/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		checkHarklorPrecursor(hardklorPath, searchPath, swathPath);
		//getSWATHPrecursors(hardklorPath, swathPath);
		
	}
	
}
