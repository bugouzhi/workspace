package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

/**
 * This class holds a series of utility static method to compute various statistics for 
 * SWATH spectra
 * @author Jian Wang
 *
 */
public class SWATHStatistics {
	//test if all peaks matched but intensity flat, distribution of projcos
	public static void testProjCosine(){
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		//Spectrum s1 = reader.getSpectrum(1631);
		System.out.println("lib size: " + lib.getSpectrumList().size());
		List<Spectrum> specList = lib.getSpectrumList();
		for(int i = 0; i < lib.getSpectrumList().size(); i++){
			Spectrum s1 = specList.get(i);
			s1.windowFilterPeaks(6, 25);
			s1.mergePeaks(s1, 0.075);
			s1.sqrtSpectrum();
			Spectrum s2 = new Spectrum(s1);
			for(int j = 0; j < s2.getPeak().size(); j++){
				s2.getPeak().get(j).setIntensity(5.0); //fixed intensity
			}
			double projcos = s1.projectedCosine(s2, 0.05);
			System.out.println("matching: " + s1.spectrumName + "\t" + s1.parentMass + "\t" + s1.charge 
					+"\tnumPeaks:\t" + s1.getPeak().size()
					+ " cosine: " + projcos);

		}
	}
	
	public static void getSWATHFiltering(){
		String queryFile = "../mixture_linked/msdata/gringar/swath_development/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
		String resultFile = "../mixture_linked/swath_test_search_win15_25_filter_pepFDR.txt";
		//MZXMLReader reader = new MZXMLReader(queryFile);
		ConsensusSpectrumReader reader = new ConsensusSpectrumReader(queryFile);
		reader.numNeighbors = 0;
		reader.minInt = 0;
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		List<String> results = Utils.FileIOUtils.createListFromFile(resultFile);
		for(int i = 0; i < results.size(); i++){
			String[] tokens = results.get(i).split("\\t");
			if(tokens.length < 6  || lib.getSpectra(tokens[4] + "." + tokens[6]) == null){
				continue; //skipping non-formatted lines
			}
			String pepKey = tokens[4] + "." + tokens[6];
			System.out.println(pepKey);
			Spectrum libEntry = lib.getSpectra(pepKey).get(0);
			int swathScan = Integer.parseInt(tokens[1]);
			Spectrum swathSpec = reader.getSpectrum(swathScan);
			//System.out.println(swathSpec);
			libEntry.windowFilterPeaks2(6, 25);
			swathSpec.filterPeaksByIntensity(50);
			swathSpec.toRelIntensity();
			swathSpec.filterPeaksByIntensity(0.03);
			//swathSpec.computePeaksZScore(0.6);
			//swathSpec.windowFilterPeaks2(15, 25);
			//swathSpec.filterPeaks(700);
			//swathSpec.filterPeaksByRankScore(4);
			double[] shared = swathSpec.projectedShare(libEntry, 0.05, 5);
			System.out.println("ZPeakStat\t"+tokens[4] + "\t" +  tokens[8] + "\t" + tokens[1] + "\t" +  tokens[7] 
			         + "\t" + swathSpec.getPeak().size() +"\t"
					+ libEntry.getPeak().size() + "\t"+ shared[0] +"\t" + shared[1]);
		}
	}
	
	//we observe low intensity peaks in SWATH and in general in the 5600 seem to have low variabilty
	//i.e. fixed intensity value, want to test to what extended this is true 
	public static void testSWATHIntensityVar(){
		String spectrumFile = "../mixture_linked/msdata/gringar/swath_development/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		//String spectrumFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		//LargeSpectrumLibIterator reader = new LargeSpectrumLibIterator(spectrumFile);
		int Bins = 501;
		double[] IntRange = new double[Bins];
		long[] count = new long[Bins];
		double[] means = new double[Bins];
		double[] variances = new double[Bins];
		double minInt = 0;
		double stepInt = 1;
		for(int i = 0; i < IntRange.length; i++){
			IntRange[i] = minInt+i*stepInt;
		}
		IntRange[Bins-1] = Double.MAX_VALUE;
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			//if(true) continue;
			//s.windowFilterPeaks2(15, 25);
			s.filterPeaksByIntensity(30);
			s.computePeaksZScore(0.5);
			List<Peak> pList = s.getPeak();
			for(int i = 0; i < pList.size(); i++){
				Peak p = pList.get(i);
				int index = ArrayUtils.getIntervalIndex(p.getIntensity(), IntRange);
				if(index < IntRange.length){
					if(count[index] > 0){	
						means[index]=(count[index]/(double)(count[index]+1))*means[index] + p.getIntensity()/(count[index]+1);
					}else{
						means[index]+= p.getIntensity();
					}
					count[index]+=1;
				}
			}
			if(s.scanNumber > 70000){
				break;
			}
		}
		reader = new MZXMLReader(spectrumFile);
		//reader = new LargeSpectrumLibIterator(spectrumFile);
		//computing var
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			List<Peak> pList = s.getPeak();
			for(int i = 0; i < pList.size(); i++){
				Peak p = pList.get(i);
				int index = ArrayUtils.getIntervalIndex(p.getIntensity(), IntRange);
				if(index < IntRange.length){
					double dev = p.getIntensity() - means[index];
					variances[index] += dev*dev / count[index];
				}
			}
			if(s.scanNumber > 70000){
				break;
			}

		}
		
		for(int i = 0; i < means.length; i++){
			System.out.println(minInt+i*stepInt + "\t" + means[i] + "\t" + variances[i] +"\t" + count[i]);
		}
		
	}
	
	//checking the varaition of masses
	public static void testSWATHMassVar(){
		String spectrumFile = "../mixture_linked/msdata/gringar/swath_development/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		//String spectrumFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		//LargeSpectrumLibIterator reader = new LargeSpectrumLibIterator(spectrumFile);
		double maxMass = 2000;
		double resolution = 0.05;
		int Bins = 40001;
		double[] MassRange = new double[Bins];
		long[] count = new long[Bins];
		double[] means = new double[Bins];
		double[] variances = new double[Bins];
		double minInt = 0;
		double stepInt = resolution;
		for(int i = 0; i < MassRange.length; i++){
			MassRange[i] = minInt+i*stepInt;
		}
		MassRange[Bins-1] = Double.MAX_VALUE;
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			s.windowFilterPeaks2(15, 25);
			//s.filterPeaksByIntensity(30);
			//s.computePeaksZScore(0.5);
			List<Peak> pList = s.getPeak();
			for(int i = 0; i < pList.size(); i++){
				Peak p = pList.get(i);
				double pmass = p.getMass()*0.9995;
				int index = (int)(pmass / resolution);
				if(pmass > maxMass){
					index = MassRange.length-1;
				}
				if(index < MassRange.length){
					if(count[index] > 0){	
						means[index]=(count[index]/(double)(count[index]+1))*means[index] + p.getIntensity()/(count[index]+1);
					}else{
						means[index]+= pmass;
					}
					count[index]+=1;
				}
			}
			if(s.scanNumber > 70000){
				break;
			}
		}
		reader = new MZXMLReader(spectrumFile);
		//reader = new LargeSpectrumLibIterator(spectrumFile);
		//computing var
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			List<Peak> pList = s.getPeak();
			for(int i = 0; i < pList.size(); i++){
				Peak p = pList.get(i);
				double pmass = p.getMass()*0.9995;
				int index = (int)(pmass / resolution);
				if(index < MassRange.length){
					double dev = p.getIntensity() - means[index];
					variances[index] += dev*dev / count[index];
				}
			}
			if(s.scanNumber > 70000){
				break;
			}

		}
		
		for(int i = 0; i < means.length; i++){
			System.out.println(minInt+i*stepInt + "\t" + means[i] + "\t" + variances[i] +"\t" + count[i]);
		}
		
	}
	
	public static void getLibraryStat(){
		String libraryFile = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
		SpectrumLib lib = new SpectrumLib(libraryFile, "MGF");
		lib.windowFilterPeaks(6, 25);
		lib.normIntensity();
		List<Spectrum> specList = lib.getSpectrumList();
		List<double[]> scores = new ArrayList<double[]>();
		for(int rank = 1; rank < 100; rank ++){
			double[] projCosines = new double[specList.size()];
			for(int i = 0; i < specList.size(); i++){
				Spectrum entry = specList.get(i);
				entry.mergePeaks(entry, 0.06);
				Spectrum copy = new Spectrum(entry);
				copy.filterPeaks(rank);
				System.out.println("number peaks: " + copy.getPeak().size());
				projCosines[i] = entry.projectedCosine(copy, 0.05);
				if(projCosines[i] > 0.85){
					//System.out.println(entry);
				}
				System.out.println("rank: " + rank + "\tprojcos: " + projCosines[i]);
				
			}
		}
	}
	
	
	
	
	
	public static void main(String[] args){
		//simpleCosTest();
		//testSWATHIntensityVar();
		//getSWATHFiltering();
		//testProjCosine();
		//testSWATHMassVar();
		getLibraryStat();
		
	}

}
