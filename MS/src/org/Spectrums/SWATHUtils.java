package org.Spectrums;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.systemsbiology.jrap.stax.Scan;

/**
 * Some utility function for SWATH data
 * @author Jian Wang
 *
 */
public class SWATHUtils {
	public static final int DALTON = 1;
	public static final int PPM = 2;
	
	
	public static boolean checkPrecursor(double ms2, double precursor, double windowWidth){
		 double diff = precursor - ms2 + 4.0;
		 return diff > 0 && diff < windowWidth;
	}
	
	public static boolean checkMass(double ms2, double precursor, double tolerance, int mode){
		return SpectrumUtil.checkMass(ms2, precursor, tolerance, mode);
	}
	
	
	public static int[] peakIntDistr(Spectrum s, double min, double max, double binWidth){
		List<Peak> sorted = new ArrayList();
		sorted.addAll(s.getPeak());
		Collections.sort(sorted, PeakIntensityComparator.comparator);
		int bins = (int)Math.ceil((max - min)/binWidth);
		int[] counts = new int[bins];
		for(int i = 0; i < sorted.size(); i++){
			
		}
		return null;
	}
	
	public static double getRT(Scan s){
		String rt = s.getHeader().getRetentionTime();
		return Double.parseDouble(rt.substring(2, rt.length()-1));
	}
	
	public static int getSWATHMS1Scan(Spectrum s){
		return (int)(s.scanNumber / 35.0)*35+1;
	}
	
	//convert a spectral library to the peakview ion library text format
	public static void getIonLibrary(String spectrumFile){
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		List<Spectrum> specList = lib.getAllSpectrums();
		int counter = 0;
		int ProteinCounter = 1; //not sure if this is protein group number but in ionlibrary seems unique to each proteins
		String header = "Q1\tQ3\tRT_detected\tprotein_name\tisotype\trelative_intensity\tstripped_sequence\tmodification_sequence\tprec_z\tfrg_type\tfrg_z\tfrg_nr\tiRT\tuniprot_id\tscore\tdecoy\tprec_y\tconfidence\tshared\tN\tmods\tnterm\tcterm";
		System.out.println(header);
		Map<String, Integer> proteinMap = new HashMap<String, Integer>();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.peptide.contains("+")){
				continue;
			}
			SimpleMatchingGraph g = SpectrumUtil.getMatchGraph(s, 0.05);
			int matchCount = 0;
			StringBuffer outbuffer = new StringBuffer();
			for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
				Peak p = it.next();
				List neighs = g.getNeighbors(p);
				LabelledPeak closestP = null;
				double smallestErr = 10000;
				for(int n = 0; n < neighs.size(); n++){
					LabelledPeak currPeak = (LabelledPeak)neighs.get(n);
					double currDiff = Math.abs(currPeak.getMass()-p.getMass());
					closestP = currDiff < smallestErr ? currPeak : closestP;
					smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
				}
				String[] tokens = s.spectrumName.split("\\s+");
				int protInd = 1;
				s.protein = tokens[7];
				if(proteinMap.containsKey(s.protein)){
					protInd = proteinMap.get(s.protein);
				}else{
					protInd = ProteinCounter++;
					proteinMap.put(s.protein, protInd);
				}
				s.rt = Double.parseDouble(tokens[5].substring(2, tokens[5].length()-1));
				if(closestP != null 
						&& (closestP.getType().equals("b") || closestP.getType().equals("y"))){
					matchCount++;
				}
				if(neighs.size() > 0){
					outbuffer.append(s.parentMass + "\t");
					outbuffer.append(p.getMass() + "\t");
					outbuffer.append(s.rt/60 + "\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(" " + "\t");
					outbuffer.append(p.getIntensity() + "\t");
					outbuffer.append(s.peptide.replaceAll("[0-9\\.\\+\\-]", "") +"\t");
					outbuffer.append(s.peptide.replaceAll("[0-9\\.\\+\\-]", "") +"\t");
					//System.out.print(s.peptide + "\t");
					outbuffer.append(s.charge  + "\t");
					outbuffer.append(closestP.getType() + "\t");
					outbuffer.append(closestP.getCharge() + "\t");
					outbuffer.append(closestP.getPos() + "\t");
					outbuffer.append(s.rt/60 + "\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(-1 + "\t");
					outbuffer.append("False" + "\t");
					outbuffer.append("10000" + "\t");
					outbuffer.append("1" + "\t");
					outbuffer.append("False" + "\t");
					outbuffer.append(protInd + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("\n");
				}
			}
			if(matchCount > 6){
				System.out.print(outbuffer);
			}else{
				//System.out.println("warning less than 10 peaks in libray spectrum");
			}
			counter++;
		}
	}
	
	public static void testGenerateIonLibrary(){
		String specLibFile = "../mixture_linked/UPS_EcoliREP3_newstock2lib.mgf";
		getIonLibrary(specLibFile);
	}
	public static void main(String[] args){
		testGenerateIonLibrary();
	}
}

