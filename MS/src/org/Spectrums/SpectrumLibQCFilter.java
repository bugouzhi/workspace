package org.Spectrums;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import Utils.SpectrumFilter;


/**
 * A set of filter designed to ensure the quality of the spectral library
 * This default set of filters are adapted from those used in NIST
 * User can add or remove filters from the default sets
 * @author Jian
 *
 */
public class SpectrumLibQCFilter {
	//instrument resolution
	public static int HI_RES = 0;
	public static int LOW_RES = 1;
	private List<SpectrumQualityFilter> filters;
	
	public SpectrumLibQCFilter(){
		this.filters = createFilterNIST();
	}
	
	public Set<Spectrum> filterLibrary(SpectrumLib lib){
		Set<Spectrum> removedSpects = new HashSet<Spectrum>();
		for(Iterator<String> it = lib.getSpectrumLibrary().keySet().iterator(); it.hasNext();){
			String key = it.next();
			List<Spectrum> specs = lib.getSpectra(key);
			for(Iterator<Spectrum> specIt = specs.iterator(); specIt.hasNext();){
				Spectrum s = specIt.next();
				for(int j = 0 ; j < filters.size(); j++){
					if(!filters.get(j).accept(s)){
						removedSpects.add(s);
						specIt.remove();
						break;
					}
				}
			}
		}
		return removedSpects;
	}
	
	/**
	 * We create the NIST standard filters for library construction
	 * @return
	 */
	public List<SpectrumQualityFilter> createFilterNIST(){
		List<SpectrumQualityFilter> filters = new ArrayList<SpectrumQualityFilter>();
		filters.add(new PepIonSignificanceFilter());
		filters.add(new PrecursorMassErrorFilter());
		filters.add(new UnIdentFragmentIonsFilter());
		filters.add(new SufficientIonsFilter());
		filters.add(new PrincipleChargeStateFilter());
		return filters;
	}
	
	
	public static class PepIonSignificanceFilter implements SpectrumQualityFilter{
		//the actual filter from NIST seems a bit complex and require external software
		//and does not seems do much more than checking for precursor presence, so here
		//we do a simple lookup of whether precursor exists
		private double tolerance = 0.05;
		//PrecursorMassChecker checker;
		
		//public PepIonSignificanceFilter(){
		//	this.checker = new PrecursorMassChecker();
		//}
		
		@Override
		public boolean accept(Spectrum s) {
			// TODO Auto-generated method stub
			//Peptide p = new Peptide(s.peptide, s.charge);
			//int matchedPrecursorPeaks = checker.matchPrecursorProfile(s.scanNumber, p, this.tolerance);
			//return matchedPrecursorPeaks > 0;
			return s.getPrecursorInt() > 0;
		}
		
		
	}
	
	public static class PrecursorMassErrorFilter implements SpectrumQualityFilter{
		private double tolerance_hiRes = 5;
		private double tolerance_lowRes = 0.25;
		private int toleranceMode;
		int resMode;
		
		@Override
		public boolean accept(Spectrum s) {
			Peptide p = new Peptide(s.peptide, s.charge);
			if(resMode == SpectrumLibQCFilter.HI_RES){
				return Mass.checkMass(s.parentMass, s.charge, this.tolerance_hiRes, Mass.DIFF_PPM);
			}else if(resMode == SpectrumLibQCFilter.LOW_RES){
				return Mass.checkMass(s.parentMass, s.charge, this.tolerance_lowRes, Mass.DIFF_DA);
			}
			return false;
		}
		
	}
	
	
	public static class UnIdentFragmentIonsFilter implements SpectrumQualityFilter{
		private double minExplInt = 1-0.32;
		private double minExpPeaks = 1-0.36;
		private double tolerance = 0.5;
		@Override
		public boolean accept(Spectrum s) {
			Spectrum copy = new Spectrum(s);
			Peptide p = new Peptide(s.peptide, s.charge);
			TheoreticalSpectrum t = new TheoreticalSpectrum(p);
			double fullExpInt = t.explainedPeaks2(t, s, this.tolerance);
			double matchedPeaksFull = s.sharePeaks(t,tolerance) / s.getPeak().size();
			copy.filterPeaks(20);
			double topExpInt = t.explainedPeaks2(t, copy, this.tolerance);
			double matchedPeakTop = s.sharePeaks(t, this.tolerance) / copy.getPeak().size();
			double subFilter1 =  (fullExpInt + topExpInt)/2;
			double subFilter2 = (matchedPeaksFull + matchedPeakTop)/2;
			return subFilter1 < this.minExplInt & subFilter2 < this.minExpPeaks;
		}
		
	}
	
	
	public static class SufficientIonsFilter implements SpectrumQualityFilter{
		private double minLargeIonFract2 = 0.2;
		private double minLargeIonFract3 = 0.3;
		private double minLargeIOnFract4plus = 0.36;
		@Override
		public boolean accept(Spectrum s) {
			int largeIonCount = 0;
			Spectrum copy = new Spectrum(s);
			copy.filterPeaks(20);
			for(int i = 0; i < copy.getPeak().size(); i++){
				Peak p = copy.getPeak().get(i);
				if(p.getMass() > s.parentMass){
					largeIonCount++;
				}
			}
			double fractLargeIons = largeIonCount/ copy.getPeak().size();
			if(s.charge == 2){
				return fractLargeIons > this.minLargeIonFract2;
			}else if(s.charge == 3){
				return fractLargeIons > this.minLargeIonFract3; 
			}else if(s.charge > 3){
				return fractLargeIons > this.minLargeIOnFract4plus;
			}
			return false;
		}
		
	}
	
	public static class PrincipleChargeStateFilter implements SpectrumQualityFilter{
		@Override
		public boolean accept(Spectrum s) {
			int chargedResCount = 0;
			String pep = s.peptide;
			for(int i = 0; i < pep.length(); i++){
				if(pep.charAt(i) == 'K' || pep.charAt(i) == 'R' || pep.charAt(i) == 'H'){
					chargedResCount++;
				}
			}
			chargedResCount++; // for nterm amine group
			return s.charge == chargedResCount;
		}
		
	}
	
	
	public static void testSpecLibQCFilter(){
		String specDir = "f:/workspace/msdata/UPS_plus_background/UPS_Ecoli/IDA_combine/";
		String psmFile = "../mixture_linked/LibraryCreation/ACG_swathdevelopment_IDA_combined_1PSMFDR_msgfdb_sorted.txt";
		String libOutFile = "../mixture_linked/LibraryCreation/testLib.mgf";
		String logOutFile = "../mixture_linked/LibraryCreation/log.txt";
		SpectrumLibConstructor constructor = new SpectrumLibConstructor(specDir, psmFile, ResultColumnIndex.MSGFDB_INDEX);
		SpectrumFilter filter = new Utils.SpectrumFilter();
		SpectrumLib lib = constructor.constructSpectralLibrary(libOutFile, logOutFile, filter, 0.05);
		SpectrumLibQCFilter qc = new SpectrumLibQCFilter();
		Set<Spectrum> removed = qc.filterLibrary(lib);
		System.out.println("Library entries did not pass filter " + removed.size());
	}
	
	public static void main(String[] args){
		testSpecLibQCFilter();
	}

}
