package org.Spectrums;

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
	int maxScan = 500000;
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
			if(s.scanNumber < minScan || s.scanNumber > maxScan){
				continue;
			}
			double parentMass = s.parentMass*s.charge-s.charge*Mass.PROTON_MASS;
			spectrumCount[ArrayUtils.getIntervalIndex(parentMass, massIntervals)]++;
			if(ArrayUtils.getIntervalIndex(parentMass, massIntervals)==2){
				System.out.println(s.scanNumber + "\t" + s.parentMass + "\t" + s.charge);
			}
		}
		this.reader = new MZXMLReader(this.msFile);
		for(int i = 1; i <= spectrumCount.length; i++){
			System.out.println("MS/MS @ mass " + massIntervals[i-1] + "--" + massIntervals[i] + "\t:\t" + spectrumCount[i-1]);
		}
	}
	
	
	public static void main(String[] args){
		String filename = "../mixture_linked/linked_peptide_library/toni/110617_Crosslink/tk110610_Nuno_start_peptide_1.mzXML";
		MS2Summary ms2sum = new MS2Summary(filename);
		ms2sum.minScan = 0;
		ms2sum.maxScan = 500000;
		ms2sum.getPrecursorChargeSummary();
		ms2sum.getParentMassSummary(new double[]{0, 1670, 2032, 5000});
//		ms2sum.getParentMassSummary(new double[]{0, 1734, 2089, 5000});
	}
}
