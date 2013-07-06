package org.Spectrums;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import Utils.FileIOUtils;

/**
 * Take in a SWATH file and a output file from Harklor
 * http://proteome.gs.washington.edu/software/hardklor/
 * Generate a IDA-like mgf file that associate all precursor
 * in the SWATH window with each SWATH MS/MS spectrum
 * @author Jian Wang
 *
 */
public class SWATHMGFGenerator {
	private String SWATHPath;
	private String HarklorPath;
	public int minCharge = 2;
	public int chargeInd = 2;
	public int mzInd = 4;
	
	public SWATHMGFGenerator(String SWATHPath, String HarklorPath){
		this.SWATHPath = SWATHPath;
		this.HarklorPath = HarklorPath;
	}
	
	private Map<Integer, List<double[]>> parseHarklorResults(){
		Map<Integer, List<double[]>> harklorResults = new HashMap();
		List<String> results = FileIOUtils.createListFromFile(this.HarklorPath);
		int currentMS1 = 0;
		List<double[]> currPrecursorList =  null;
		for(int i = 0; i < results.size(); i++){
			String[] tokens = results.get(i).split("\\s+");
			//System.out.println("line is : " + results.get(i));
			if(tokens[0].equals("S")){
				if(currPrecursorList != null){
					//System.out.println("adding scan: " + currentMS1 + "\t" + currPrecursorList.size());
					harklorResults.put(currentMS1, currPrecursorList);
				}
				currentMS1 = Integer.parseInt(tokens[1]);
				currPrecursorList = new ArrayList();
			}
			if(tokens[0].equals("P")){
				double mz = Double.parseDouble(tokens[mzInd]);
				double charge = Double.parseDouble(tokens[chargeInd]);
				currPrecursorList.add(new double[]{mz, charge});
			}
		}
		return harklorResults;
	}
	
	public void generateExpandedMGF(String outFile){
		Map<Integer, List<double[]>> harklorResults = parseHarklorResults();
		MZXMLReader reader = new MZXMLReader(this.SWATHPath);
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			while(reader.hasNext()){
				Spectrum s = reader.next();
				//s.filterPeaks(500);
				s.windowFilterPeaks2(15, 25);
				Spectrum copy = new Spectrum(s);
				//System.out.println("Number of peaks: " + s.getPeak().size());
				int MS1 = reader.getPrevScan(s.scanNumber, 1);
				System.out.println("Processing: " + MS1);
				if(harklorResults.containsKey(MS1)){
					List<double[]> precursorList = harklorResults.get(MS1);
					//System.out.println("MS1 is: " + MS1 + "\t" + precursorList.size());
					for(int i = 0; i < precursorList.size(); i++){
						double[] curr = precursorList.get(i);
						//System.out.println(MS1 +":\t" + s.scanNumber + "\t" + s.parentMass + "\t" + curr[0]);
						if(checkPrecursor(s.parentMass, curr[0], 25.0)
								&& curr[1] >= this.minCharge){
							//System.out.println("founded");
							copy.parentMass = curr[0];
							copy.charge = (int)curr[1];
							out.write(copy.toString());
						}
					}
				}
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	private boolean checkPrecursor(double ms2, double precursor, double windowWidth){
		 double diff = precursor - ms2 - 4.0;
		 return diff > 0 && diff < windowWidth;
	}
	
	public static void testGetMGF(){
		String harklorFile = "..//mixture_linked/Hardklor/win32/14344_UPS1_400fm_Ecolilysate_SWATH_5600_default_0.7min.hk";
		String swathFile = "..//mixture_linked/msdata/gringar/swath_development/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		SWATHMGFGenerator gen = new SWATHMGFGenerator(swathFile, harklorFile);
		gen.generateExpandedMGF("..//mixture_linked//swath_expanded.mgf");
	}
	
	public static void main(String[] args){
		testGetMGF();
	}
	
}
