package org.Spectrums;

import java.util.Iterator;

/**
 * We test a hybrid method that utilize both features from spectral matches and 
 * DBsearch to discriminate true matches from false ones
 * @author Jian Wang
 *
 */
public class MSPLITHybrid extends SpectrumLib{
	private static final long serialVersionUID = 147324738947L;
	public static void hybridSearch(SpectrumComparator comp1, SpectrumComparator comp2, Iterator<Spectrum> query, SpectrumLib lib, 
			double parentTolerance, double fragmentTolerance){
		lib.toNormVector(fragmentTolerance*2, fragmentTolerance,  2000);
		lib.normIntensity();
		lib.setParentMassTolerance(parentTolerance);
		for(Iterator<Spectrum> mixIter = query; mixIter.hasNext();){
			Spectrum mixture = mixIter.next();
			mixture.windowFilterPeaks(10, 25);
			//mixture.scaleMass(0.9995);
			Spectrum mixture2 = mixture.toNormVector(fragmentTolerance*2, fragmentTolerance, 2000);
			mixture2.sqrtSpectrum();
			mixture2 = mixture2.toNormVector(fragmentTolerance*2, fragmentTolerance, 2000);
			Spectrum best = lib.searchAndBoundLib(mixture2);
			String[] peps = best.peptide.split(" & ");
			mixture.scaleMass(1.0005);
			mixture.computePeakRank();
			if(peps[0].equals("X")){
				peps[0] = peps[0].substring(2);
			}
			if(peps[0].equals("X")){
				peps[1] = peps[1].substring(2);
			}
			SpectrumUtil.getMixtureMatchStat(mixture, peps[0], peps[1], comp1, comp2);
		}
	}
	
	public static void testHybridSearch(){
		String filename = "..\\MSPLib\\Lib\\human.msp";
		String training1 = "..\\mixture_linked\\mixtures100000_alpha0.5.mgf";
		String training2 = "..\\MSPLib\\Lib\\yeast.msp";
		String libFile = "..\\mixture_linked\\yeast_data\\yeast_decoy.sptxt";
		String mixFile = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mgf";
        SpectrumComparator comp1 = SpectrumUtil.getMixtureScorer(training1);
      //  SpectrumLib lib1 = new SpectrumLib(filename, "MSP");//
      //  lib1.removeModSpectra();
      //  lib1.scaleSpectrumMass(0.9995);
      //  SpectrumLib lib = new SpectrumLib(training2, "MSP");//
      //  lib.removeModSpectra();
      //  lib.scaleSpectrumMass(0.9995);
        SpectrumComparator comp2 = SpectrumUtil.getRankBaseScorer(training2);
        SpectrumLib lib = new SpectrumLib(libFile, "splib");
        lib.scaleSpectrumMass(0.9995);
        lib.removeModSpectra();
        //SpectrumLib lib2 = lib1.Divide();
        //SpectrumLib mixlib = lib2.createRandomMix(lib, 2000, 0.2, 0, 1, 3, false);
        //lib1=null;
        //lib2=null;
        //lib = null;
        LargeSpectrumLibIterator iter = new LargeSpectrumLibIterator(mixFile); 
        //hybridSearch(comp1, comp2, mixlib.getAllSpectrums().iterator(), lib, 3, 0.5);
        hybridSearch(comp1, comp2, iter, lib, 3, 0.5);
	}
	
	public static void main(String[] args){
		testHybridSearch();
	}
}
