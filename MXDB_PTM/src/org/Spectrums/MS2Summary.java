package org.Spectrums;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * print out some summary for mzxml file
 * @author Jian Wang
 *
 */
public class MS2Summary {
	private String msFile;
	MZXMLReader reader;
	int minScan = 0;
	int maxScan = 5000000;
	public MS2Summary(String file){
		this.msFile = file;
		reader = new MZXMLReader(file);	
	}
	
	public void getPrecursorChargeSummary(){
		Map<Integer, Integer> chargeCount = new HashMap();
		for(;reader.hasNext();){
			Spectrum s = reader.next();
			if(s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}
			if(chargeCount.containsKey(s.charge)){
				int count = chargeCount.get(s.charge);
				count++;
				chargeCount.put(s.charge, count);
			}else{
				chargeCount.put(s.charge, 1);
			}
		}
		this.reader = new MZXMLReader(this.msFile);
		for(Iterator<Integer> it = chargeCount.keySet().iterator(); it.hasNext();){
			int charge = it.next();
			System.out.println("MS/MS @ charge " + charge + " :\t" + chargeCount.get(charge));
		}
	}
	
	public void getParentMassSummary(double[] massIntervals){
		int[] spectrumCount = new int[massIntervals.length-1];
		for(;reader.hasNext();){
			Spectrum s = reader.next();
			if(s.scanNumber < minScan || s.scanNumber > maxScan || s.charge < 0){
				continue;
			}
			double parentMass = s.parentMass*s.charge-s.charge*Mass.PROTON_MASS;
			spectrumCount[ArrayUtils.getIntervalIndex(parentMass, massIntervals)]++;
			if(ArrayUtils.getIntervalIndex(parentMass, massIntervals)==2){
				//System.out.println(s.scanNumber + "\t" + s.parentMass + "\t" + s.charge);
			}
		}
		this.reader = new MZXMLReader(this.msFile);
		for(int i = 1; i <= spectrumCount.length; i++){
			System.out.println("MS/MS @ mass " + massIntervals[i-1] + "--" + massIntervals[i] + "\t:\t" + spectrumCount[i-1]);
		}
	}
	
	public static void getSummaryForMSFiles(String directory){
		File dir = new File(directory);
		File[] files = dir.listFiles();
		for(int i =0; i < files.length; i++){
			File current = files[i];
			//System.out.println("file is: " + current.getName());
			String name = current.getName();
			if(name.matches(".+\\.mzXML")){
				System.out.println("MS data file:\t" + name);
				MS2Summary ms2sum = new MS2Summary(directory+"/"+name);
				ms2sum.minScan = 0;
				ms2sum.maxScan = 15000000;
				ms2sum.getPrecursorChargeSummary();
				ms2sum.getParentMassSummary(new double[]{0, 1187, 1788, 2375, 15000});
			}
		}
	}
	
	public static void main(String[] args){
		String filename = "../mixture_linked/msdata/linked_peptide_library/ACG_disulfide_library/Orbi_Elite/PepLib1_300ng_trp_Elite_CID_25rep.mzXML";
		String dir = "../mixture_linked/msdata/toni_110510_Crosslink/";
//		getSummaryForMSFiles(dir);
		MS2Summary ms2sum = new MS2Summary(filename);
//		ms2sum.minScan = 0;
//		ms2sum.maxScan = 15000000;
		ms2sum.getPrecursorChargeSummary();
//		ms2sum.getParentMassSummary(new double[]{0, 1661, 2318, 3322, 15000});  //library_crosslink 9 undigested
		//ms2sum.getParentMassSummary(new double[]{0, 1402, 2057, 2800, 15000});
//		ms2sum.getParentMassSummary(new double[]{0, 1446, 2049, 2800, 15000}); //library_crosslink 9  digested
//		ms2sum.getParentMassSummary(new double[]{0, 1187, 1788, 2375, 15000});
		//ms2sum.getParentMassSummary(new double[]{0, 1735, 2239, 15000});
		ms2sum.getParentMassSummary(new double[]{0, 1735-300, 2239, 2700, 3500, 15000});
		
	}
}
