package org.Spectrums;
/**
 * Do some preprocessing for SWATH spectrum
 * it seems the centroided SWATH spectrum is not very good
 * some information may be lost in its native peak-picking form
 * @author Jian
 *
 */
public class SWATHPreprocessor {
	public static void testCompareRawToPeakPicked(){
		String raw = "f:/workspace/msdata/APSWATH_EIF4A2_MEPCE_GFP_set/Swath_EIF4aJune7-Biorep2.mzXML";
		String peakPicked = "";
		MZXMLReader reader = new MZXMLReader(raw);
		int startscan = 20000;
		for(int i = startscan; i < 30000; i++){
			Spectrum s = reader.getSpectrum(i);
			System.out.println("Full spectrum with peaks: " + s.getPeak().size());
			s.mergePeaks(s, 0.03);
			System.out.println("Merged spectrum with peaks: " + s.getPeak().size());
			s.windowFilterPeaks2(15, 25);
			System.out.println("Filtered spectrum with peaks: " + s.getPeak().size());
		}
		
	}
	
	
	public static void main(String[] args){
		testCompareRawToPeakPicked();
	}
	

}
