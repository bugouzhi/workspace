package org.Spectrums;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;

import IO.MZXMLReader;
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
	
	public void generateExpandedMGF(String outFile){
		HarklorParser parser = new HarklorParser(this.HarklorPath);
		Map<Integer, SortedMap<Double, Double[]>> harklorResults = parser.precursorList;
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
					SortedMap<Double, Double[]> precursorList = harklorResults.get(MS1);
					//System.out.println("MS1 is: " + MS1 + "\t" + precursorList.size());
					Iterator<Double[]> precursors = precursorList.subMap(s.parentMass -5, s.parentMass + 26).values().iterator();
					while(precursors.hasNext()){
							Double[] curr = precursors.next();
							copy.parentMass = curr[0];
							copy.charge = (int)(curr[1].doubleValue());
							out.write(copy.toString());
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
