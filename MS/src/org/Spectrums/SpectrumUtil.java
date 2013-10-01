package org.Spectrums;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import Utils.FileIOUtils;



/**
 * Various utility method for spectrum 
 * @author jian wang
 *
 */
public class SpectrumUtil {
	public static final int DALTON = 1;
	public static final int PPM = 2;
	
	public static Spectrum toVector(Spectrum s){
		return null;
	}
	
	public static double magnitude(Spectrum s){
		return 0.0;
	}
	
	public static double sumOfPeaks(Spectrum s){
		return 0.0;
	}
	
	public static void normalize(Spectrum s){
		
	}
	
	public static void normalizeByTotalIntensity(Spectrum s){
		
	}
	
	public static void scaleSpectrum(double scale){
		
	}
	
	public static void shiftSpectrum(double shift){
		
	}
	
	public static void toRelIntensity(){
		
	}
	
	public static void sqrtSpectrum(){
		
	}
	
	
	public static boolean checkMass(double ms2, double precursor, double tolerance, int mode){
		 double diff = precursor - ms2;
		 diff = Math.abs(diff);
		 
		//System.out.println("diff: " + diff + "\t" + diff*1000000/precursor +"\t" + tolerance + "\t" + (diff < tolerance));

		 if(mode == DALTON){
			 return diff < tolerance;
		 }
		 if(mode == PPM){
			 return diff*1000000/precursor < tolerance;
		 }
		 return false;
	}
	
	public static List<Peak> copyPeaks(Spectrum s){
		List<Peak> copies = new ArrayList<Peak>();
		for(Iterator<Peak> it = s.getPeaks().iterator(); it.hasNext();){
			Peak current = it.next();
			Peak copy = new Peak(current);
			copies.add(copy);
		}
		return copies;
	}
	
	public static List<Peak> scaleIntensity(List<Peak> pList, double scale){
		for(Iterator<Peak> it = pList.iterator(); it.hasNext();){
			Peak current = it.next();
			current.setIntensity(current.getIntensity()*scale);
		}
		return pList;
	}
	
	public static List<Peak> toRelativeIntensity(List<Peak> pList){
		double maxIntensity=0;
		for(int i = 0; i < pList.size(); i++){
			double intensity = pList.get(i).getIntensity();
			maxIntensity = maxIntensity > intensity ? maxIntensity : intensity;
		}
		
		return scaleIntensity(pList, 1/maxIntensity);
	}
	//normalized peak rank, range from 1 to 100
	public static int normalizedRank(Peak p, Spectrum s){
		return normalizedRank(p, s.getPeak().size());
	}
	
	public static int normalizedRank(Peak p, int max){
		double relrank = p.getRank()/(double)max*100.0;
		return (int)Math.round(relrank);
	}
	
	public static List<Spectrum> getRandomSpectrum(SpectrumLib lib, int size){
		List<Spectrum> specList = new ArrayList();
		Map<Spectrum, Spectrum> map = new HashMap();
		int count = 0;
		while(count < size){
			Spectrum s = lib.getRandomSpectrum();
			if(!map.containsKey(s)){
				specList.add(s);
				count++;
			}
		}
		return specList;
	}
	
	public static List<Spectrum> getRandomSpectrum(SpectrumLib lib, int size, int charge){
		List<Spectrum> specList = new ArrayList();
		Map<Spectrum, Spectrum> map = new HashMap();
		int count = 0;
		while(count < size){
			Spectrum s = lib.getRandomSpectrum();
			if(!map.containsKey(s) && s.charge <= charge){
				specList.add(s);
				count++;
			}
		}
		return specList;
	}
	
	public static List<Spectrum> getSpectra(String file, SpectrumLib lib){
		Map<String, String> spectrumIds = 
			FileIOUtils.createTableFromFile(file, 0, 0);
		List<Spectrum> list = new ArrayList();
		for(Iterator<String> idIter = spectrumIds.keySet().iterator(); idIter.hasNext();){
			String key = idIter.next();
			if(lib.getSpectrumLibrary().containsKey(key)){
				Spectrum s = lib.getSpectra(key).get(0);
				list.add(s);
			}
		}
		return list;
	}
	
	public static void computePeakRanks(List<Peak> peaks){
		List<Peak> sortedList = new Vector<Peak>();
		sortedList.addAll(peaks);
		Collections.sort(sortedList, PeakIntensityComparator.comparator);
		Map<Peak, Peak> peakMap = new HashMap<Peak, Peak>();
		Peak p;
		for(int i = 0, size = sortedList.size(); i < size; i++){
			p = sortedList.get(i);
			p.setRank(size - i);
		}
	}
	
	public static void printLibPeptides(){
		String file = ".\\MSPLib\\Lib\\yeast.msp";
		SpectrumLib lib = new SpectrumLib(file, "MSP");
		lib.removeModSpectra();
		List<Spectrum> specList = lib.getSpectrumList();
		for(int i = 0; i < specList.size(); i++){
			System.out.println(specList.get(i).getPeptide());
		}
		
	}
	
	public static void analyzeMixtureSpectralMatches(){
		String spectrumFile = ".\\MSPLib\\Lib\\yeast.msp";
		String mixSpectrumFile = "..\\mixture_linked\\yeast_data\\klc_122007p_yeast_digest1.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		lib1.removeModSpectra();
		lib1.scaleSpectrumMass(0.9995);
		System.out.println("Done loading spectrum library");
		SpectrumLib mix = new SpectrumLib(mixSpectrumFile, "MGF");
		System.out.println("Done loading mixture library");
		String tripletFile = "..\\mixture_linked\\mixtureMatches.txt";
		mix.scaleSpectrumMass(0.9995);
		mix.toRelativeIntensity();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); 
			String[] tokens, tokens2;
			Spectrum s1, s2, m;
			while(currentLine != null){
				tokens = currentLine.split("\t");
				m = mix.getSpectra(tokens[0]).get(0);
				m.windowFilterPeaks(10, 25);
				m = m.toNormVector();
				m.sqrtSpectrum();
				s1 = lib1.getSpectra(tokens[1]).get(0);
				//s1.scaleMass(0.9995);
				s1 = s1.toNormVector();
				s1.sqrtSpectrum();
				s2 = lib1.getSpectra(tokens[2]).get(0);
				//s2.scaleMass(0.9995);
				s2 = s2.toNormVector();
				s2.sqrtSpectrum();
				double alpha = m.alpha(s1, s2);
				System.out.println("Mathicng: " + m.peptide + " to " + s1.peptide + " and " + s2.peptide + " has similarity: " + m.maxScore(s1, s2, alpha)
						+ "\t" + m.cosineSim(s1) + "\t" + m.cosineSim(s2) + "\t" + s1.cosineSim(s2) + "\t" + alpha + "\t" + m.residual(s1));
				currentLine = bf.readLine();
			}
			
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void getMixtureMatchStat(Spectrum query, String p1, String p2, 
		SpectrumComparator comp1, SpectrumComparator comp2){
		TheoreticalSpectrum t1 = new TheoreticalSpectrum(p1);
		TheoreticalSpectrum t2 = new TheoreticalSpectrum(p2);
		double score1 = comp2.compare(t1, query);
		double score2 = comp2.compare(t2, query);
		TheoreticalSpectrum th = new TheoreticalSpectrum(p1, p2);
		double cscore = comp1.compare(th, query);
		if(p1.equals(p2)){ //we do not allow same peptides 
			return;
		}
		double[] stat = th.analyzeMixtureAnnotation(query, p1, p2);
		String bestpeptide = p1 + " & " + p2;
		System.out.print("Spectrum: " + query.getPeptide() + " has best match: " +  bestpeptide 
				+ " with score:\t" +  cscore + "\t"  + score1 + "\t" + score2 + "\t" + score1/p1.length() + "\t" + score2/p2.length()
				+ "\t" + stat[0] + "\t" + stat[1] + "\t"
				+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
				+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
				+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12]);
		//System.out.println("\t" + checkPeptidepair(bestpeptide, query.peptide));	
		System.out.println();
		
	}
	
	public static SpectrumComparator getMixtureScorer(String trainingFile){
		System.out.println("Starting mixture-training");
		MixturePeakScoreLearner peakscorer3 = new MixturePeakScoreLearner(trainingFile);
		peakscorer3.getMixtureIonCount();
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(peakscorer3);
		return scorer;
	}
	

	
	public static SpectrumComparator getLPeakRankBaseScorer(String file){
		LPeakRankBaseScorer learner = null;
		if(file.contains(".msp")){
			SpectrumLib lib1 = new SpectrumLib(file, "MSP");
			lib1.removeModSpectra();
			lib1.computeRank();
			learner = new LPeakRankBaseScorer(lib1);
		}
		
		if(file.contains(".o")){
			learner = LPeakRankBaseScorer.loadComparator(file);
		}
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		scorer.includeNoise = false;
		return scorer;
	}
	
	
	
	
	public static SpectrumComparator getLPeakRankBaseScorer2(String trainingFile){
		LPeakRankBaseScorer learner = new LPeakRankBaseScorer(trainingFile);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		return scorer;
	}
	
	public static SpectrumComparator getLMixtureScorer(String trainingFile){
		LMixturePeakScoreLearner learner = null;
		if(trainingFile.contains(".o")){
			learner = LMixturePeakScoreLearner.loadComparator(trainingFile);
		}else{
			learner = new LMixturePeakScoreLearner(trainingFile);
		}
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(learner);
		return scorer;
	}
	
	public static SpectrumComparator getLinkedPeptideScorer(String file){
		LinkedPeptidePeakScoreLearner learner = null;
		if(file.contains(".mgf")){
			learner = new LinkedPeptidePeakScoreLearner(file);
			learner.getLinkedIonCount();
		}
		if(file.contains(".o")){
			learner = LinkedPeptidePeakScoreLearner.loadComparator(file);
		}
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(learner);
		return scorer;
	}
	
	public static SpectrumComparator getLinkedSUMOScorer(String file){
		LinkedPeptidePeakScoreLearner learner = null;
		if(file.contains(".mgf")){
			learner = new LinkedPeptidePeakScoreLearner(file);
			learner.getLinkedIonCount();
		}
		if(file.contains(".o")){
			learner = LinkedPeptidePeakScoreLearner.loadComparator(file);
		}
		//MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(learner);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		scorer.setMinMatchedPeak(0);
		return scorer;
	}
	
	public static SpectrumComparator getLinkedPeptideSingleScorer(String file){
		LinkedPeptidePeakScoreLearner learner = null;
		if(file.contains(".mgf")){
			learner = new LinkedPeptidePeakScoreLearner(file);
			learner.getLinkedIonCount();
		}
		if(file.contains(".o")){
			learner = LinkedPeptidePeakScoreLearner.loadComparator(file);
		}
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(learner);
		return scorer;
	}
	
	public static SpectrumComparator getLinkedPeptideScorer(SpectrumLib lib){
		LinkedPeptidePeakScoreLearner learner = new LinkedPeptidePeakScoreLearner(lib);
		learner.getLinkedIonCount();
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(learner);
		return scorer;
	}
	
	public static SpectrumComparator getRankBaseScorer(String trainingFile){
		SpectrumLib lib1 = null;
		if(trainingFile.endsWith("msp")){
			lib1 = new SpectrumLib(trainingFile, "MSP");
		}else if(trainingFile.endsWith("mgf")){
			lib1 = new SpectrumLib(trainingFile, "MGF");
		}
		return getRankBaseScorer(lib1);
	}
	
	
	public static PeakComparator getRankBasePeakComparator(String trainingFile){
		SpectrumLib lib1 = new SpectrumLib(trainingFile, "MSP");
		return getRankBasePeakComparator(lib1);

	}
	
	public static SpectrumComparator getRankBaseScorer(SpectrumLib lib){
		PeakComparator peakscorer2 = getRankBasePeakComparator(lib);
		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);
		//SimpleProbabilisticScorer2 scorer1 = new SimpleProbabilisticScorer2(peakscorer2);
		return scorer1;
	}
	
	public static PeakComparator getRankBasePeakComparator(SpectrumLib lib){
		lib.removeModSpectra();
		lib.computeRank();
		System.out.println("Starting rank-base training");
		RankBaseScoreLearner peakscorer2 = new RankBaseScoreLearner(lib);
		peakscorer2.getIonsCount();
		return peakscorer2;
	}
	
	public static SpectrumComparator getSimpleScorer(String trainingFile, boolean includeNoise){
		SpectrumLib lib1 = new SpectrumLib(trainingFile, "MSP");
		lib1.removeModSpectra();
		lib1.computeRank();
		String probFile = ".//data//IonsScore.txt";
		String noiseModel = ".//data//NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		return filter;
	}
	
	public static void convertMSPToMGF(String filename){
	    SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
	    lib1.printLibToFile("..\\mixture_linked\\yeast_specLib.mgf", lib1);
	}
	
	public static void loadAnnotationFromFile(String filename, SpectrumLib lib){
		loadAnnotationFromFile(filename, lib.getAllSpectrums());
	}
	
	public static void loadAnnotationFromFile(String filename, List<Spectrum> specList){
		Map<String, String> annotationMap = 
			FileIOUtils.createTableFromFile(filename, 0, 1);
		int count =0;
		for(Iterator<Spectrum> it = specList.iterator(); it.hasNext();){
			Spectrum s = it.next();
			if(annotationMap.containsKey(""+s.scanNumber)){
				System.out.println("scan number: " + s.scanNumber);
				String annotation = annotationMap.get(""+s.scanNumber);
				String[] peps = annotation.split(" & ");
				s.peptide = annotation;
				count++;
				//Peptide p = new Peptide(annotation);
				//LinkedPeptide p = new LinkedPeptide(peps[0]+"--"+peps[1], s.charge,  3, 6);
				//s.setPeptide(p.getPeptide());//+"."+s.charge);
				//s.charge = p.getCharge();
				//s.parentMass = p.getParentmass();
			}
			
		}
		System.out.println("annotated: " + count);
	}

	public static void removeUnannotatedSpectra(SpectrumLib lib){
		for(Iterator<Spectrum> it = lib.getAllSpectrums().iterator(); it.hasNext();){
			Spectrum s = it.next();
			if(s.peptide == null ){
				lib.removeSpectrum(s.spectrumName);
			}else if(s.peptide.contains("Scan Number")){
				lib.removeSpectrum(s.peptide);
			}
		}
	}
	
	/**
	 * Compute mass difference both in Dalton and in ppm
	 * @param mass1
	 * @param mass2
	 * @param mode - 1 is difference in dalton, 2 is ppm
	 * @return
	 */
	public static double massDiff(double mass1, double mass2, int mode){
		if(mode == 1){
			return Math.abs(mass1 - mass2);
		}else{
			return Math.abs((mass1 - mass2)*1000000 / mass2); 
		}

	}
	
	public static SimpleMatchingGraph getMatchGraph(Spectrum s, double tolerance){
		Peptide pep = new Peptide(s.peptide, s.charge);
		TheoreticalSpectrum.prefixIons = new String[]{"b"};//, "b-H20", "b-NH3"};//, "b(iso)"};
		TheoreticalSpectrum.suffixIons = new String[]{"y"};//, "y-H20", "y-NH3"};//, "y(iso)"};
		TheoreticalSpectrum th = new TheoreticalSpectrum(pep);
		//TheoreticalSpectrum.addIsotopicPeaks(th, 1);
		SimpleMatchingGraph g = th.getMatchGraph(s, tolerance);
		return g;
	}
	
	public static void annotateSpectrumLib(){
		String libFile = "..\\mixture_linked\\yeast_data/klc_122007p_yeast_digest1.mgf";
		String annotationFile = "..\\mixture_linked\\yeast1_MSPLIT_result.txt";
		SpectrumLib lib1 = new SpectrumLib(libFile, "MGF");
		loadAnnotationFromFile(annotationFile, lib1);
		lib1.printLibToFile("..\\mixture_linked\\yeast1_annotatedLib.mgf", lib1);
	}
	
	public static void annotateSpectrumLibFromMzXML(){
		String libFile = "..\\mixture_linked\\msdata/mcl1/MCL1_Bo/MCL1_Bo_Light_200fm_40%Collision.mzXML";
		String annotationFile = "..\\mixture_linked\\testAnnotation.txt";
		String outfile = "..\\mixture_linked\\test.mgf";
		MZXMLReader reader = new MZXMLReader(libFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		loadAnnotationFromFile(annotationFile, specList);
		printSpectraToFile(outfile, getAnnotatedSpectra(specList));
	}
	
	public static void annotateSpectrumLibFromMzXMLs(String spectrumDir, String annotationFile, String outfile){
		spectrumDir = "../mixture_linked//msdata/UPS12_Human/";
		annotationFile = "..\\mixture_linked/t0";
		outfile = "..\\mixture_linked\\test.mgf";
		List<Spectrum> specList = new ArrayList<Spectrum>();
		Map<String, MZXMLReader> readers = new HashMap<String, MZXMLReader>();
		try{
			BufferedReader buff = new BufferedReader(new FileReader(annotationFile));
			BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
			String line;
			String prevFile= "";
			MZXMLReader reader = null;
			int counter = 1; //one base index
			int pepInd = 4;
			int chargeInd = 6;
			while((line=buff.readLine()) != null){
				String[] tokens = line.split("\\t");
				String file = tokens[0];
				file = file.replaceAll("\\s+", "");
				file = FileIOUtils.getFileName(file);
				if(!readers.containsKey(file)){
					MZXMLReader newreader = new MZXMLReader(spectrumDir+"\\" + file);
					readers.put(file, newreader);
				}
				reader = readers.get(file);
				System.out.println("line\t"+line);
				Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[1]));
				String peptide = tokens[pepInd];
				int charge = Integer.parseInt(tokens[chargeInd]);
				//peptide = peptide.substring(2, peptide.length()-2);
				Peptide p = new Peptide(peptide, charge);
				s.peptide=peptide;
				s.charge = charge;
				//s.protein = tokens[8];
//				if(Math.abs(p.getParentmass() - s.parentMass) > 1.2){
//					s = reader.getSpectrum(Integer.parseInt(tokens[1])-1); //try to fix off-1 difference in scan number in proteowiz conversion
//				}
//				if(Math.abs(p.getParentmass() - s.parentMass) > 1.2){
//					s = reader.getSpectrum(Integer.parseInt(tokens[1])-2); //try to fix off-1 difference in scan number in proteowiz conversion
//				}
//				if(Math.abs(p.getParentmass() - s.parentMass) > 1.2){
//					s = reader.getSpectrum(Integer.parseInt(tokens[1])-6); //try to fix off-1 difference in scan number in proteowiz conversion
//				}
//				if(Math.abs(p.getParentmass() - s.parentMass) > 1.2){
//					s = reader.getSpectrum(Integer.parseInt(tokens[1])-4); //try to fix off-1 difference in scan number in proteowiz conversion
//				}
				if(Math.abs(p.getParentmass() - s.parentMass) > 1.2){
					//System.out.println(s);
					System.out.println("warning: parent mass not matching");
				}
				s.peptide=peptide;
				s.charge = Integer.parseInt(tokens[6]);
				s.scanNumber = counter++;
				s.spectrumName = s.spectrumName + " PROTEIN: " + tokens[8]; 
				out.write(s.toString());
				out.write("\n");
				if(counter > 50){
					//break;
				}
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		//printSpectraToFile(outfile, specList);
	}
	
	public static void removeAnnotateFromMzXMLs(){
		String spectrumDir = "../mixture_linked/msdata/UPS_Ecoli";
		String annotationFile = "..\\mixture_linked/t000";
		String outfile = "..\\mixture_linked\\test.mgf";
		List<Spectrum> specList = new ArrayList<Spectrum>();
		int queryInd = 0;
		int scanInd = 1;
		int pepInd1 = 4;
		int pepInd2 = 6;
		int chargeInd1 = 6;
		Map<String, MZXMLReader> readers = new HashMap<String, MZXMLReader>();
		try{
			BufferedReader buff = new BufferedReader(new FileReader(annotationFile));
			BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
			String line;
			String prevFile= "";
			MZXMLReader reader = null;
			int counter = 1; //one base index
			while((line=buff.readLine()) != null){
				String[] tokens = line.split("\\t");
				String file = tokens[queryInd];
				int index = file.lastIndexOf("/");
				file = file.substring(index);
				file = file.replaceAll("\\s+", "");
				if(!readers.containsKey(file)){
					MZXMLReader newreader = new MZXMLReader(spectrumDir+"\\" + file);
					readers.put(file, newreader);
				}
				reader = readers.get(file);
				System.out.println("line\t"+line);
				Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[scanInd]));
				//s.windowFilterPeaks2(15, 25);
				String peptide = tokens[pepInd1];
				//peptide = peptide.substring(2, peptide.length()-2);
				//int charge = Integer.parseInt(peptide.substring(peptide.length()-1));
				int charge = Integer.parseInt(tokens[chargeInd1]);
				peptide = peptide.substring(0, peptide.length()-2);
				Peptide pep = new Peptide(peptide, charge);
				//s.peptide=peptide;
				//s.charge = charge;
				TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
				t.addIsotopicPeaks(t, 1);
				SpectrumUtil.removeAnnotatedPeaks(s, t, 0.03, 50);
				if(Math.abs(pep.getParentmass() - s.parentMass) > 1.2){
					System.out.println("warning: parent mass not matching");
				}
				//s.spectrumName = s.spectrumName + " PROTEIN: " + tokens[8]; 
				out.write(s.toString());
				out.write("\n");
				if(counter > 50){
					//break;
				}
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		//printSpectraToFile(outfile, specList);
	}
	
	public static void annotateSpectrumLibFromMGF(){
		String libFile = "..\\mixture_linked\\spectral_library\\iPRG2012_nd.mgf";
		String annotationFile = "..\\mixture_linked\\MSPLIT\\iPRG2012_nd_combined_qtop_IT_MSPLIT_IDs.txt";
		String outfile = "..\\mixture_linked\\test.mgf";
		LargeSpectrumLibIterator<Spectrum> it = new LargeSpectrumLibIterator(libFile);
		Map table = FileIOUtils.createTableFromFile(annotationFile, 0, 1);
		List<Spectrum> specList = new ArrayList<Spectrum>();
		int index = 1;
		for(;it.hasNext();){
			Spectrum s = it.next();
			String key = ""+index;
			if(table.containsKey(key)){
				System.out.println(s.spectrumName + "\t" + s.scanNumber + "\t" + s.peptide);
				s.peptide = (String)table.get(key);
				System.out.println("annotated " + key + " with: " + s.peptide);
			}
			specList.add(s);
			index++;
		}
		
		//List<Spectrum> specList = it.readAllSpectra();
		//loadAnnotationFromFile(annotationFile, specList);
		printSpectraToFile(outfile, getAnnotatedSpectra(specList));
	}
	
	
	public static List<Spectrum> getAnnotatedSpectra(List<Spectrum> specList){
		List<Spectrum> annotatedList = new ArrayList<Spectrum>();
		for(Iterator<Spectrum> it = specList.iterator(); it.hasNext();){
			Spectrum curr = it.next();
			if(!(curr.peptide.contains("DUMMY") || curr.peptide.contains("Scan Number"))){
				annotatedList.add(curr);
			}
		}
		System.out.println("size of annotated: " + annotatedList.size());
		return annotatedList;
		
	}
	
	public static void printSpectraToFile(String outfile, List<Spectrum> specList){
		try{
			BufferedWriter bo = new BufferedWriter(new FileWriter(outfile));
			Iterator<Spectrum> it = specList.iterator();
			Spectrum curr;
			while(it.hasNext()){
				bo.write(it.next().toString());
				bo.write("\n");
			}
			bo.flush();
			bo.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		}
		
	}

	
	//insert parentmass info into inspect results
	public static void insertParentMassInfo(){
		String annotationFile = "..\\mixture_linked\\inspect_result_longrun.txt";
		String spectrumFile = "..\\mixture_linked\\new80min.mgf";
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		List<String> spectrumMass = 
			FileIOUtils.createListFromFile(annotationFile);
		//System.out.println("read in " + spectrumMass.size() + " lines");
		for(Iterator<String> idIter = spectrumMass.iterator(); idIter.hasNext();){
			String line = idIter.next();
			if(line.contains("#SpectrumFile")){
				line = idIter.next(); //skipping header
			}
			String[] tokens = line.split("\t");
			Spectrum s = lib.getSpectrumById("spec_"+tokens[1]+".dta.");
			System.out.println(line + "\t" + s.parentMass);
		}
	}
	
	public static SimpleMatchingGraph constructMatchingGraph(Spectrum s1, Spectrum s2, double tolerance){
			List<Peak> list1 = new ArrayList<Peak>();
			List<Peak> list2 = new ArrayList<Peak>();
			List<Peak> pList1 = s1.getPeak(); 
			List<Peak> pList2 = s2.getPeak();
			SimpleMatchingGraph matchingGraph = new SimpleMatchingGraph();
			int i = 0, j = 0, p = 0;
			double m1, m2, diff;
			Peak p1, p2;
			for(int c = 0; c < pList2.size();c++){
				matchingGraph.addVertex(pList2.get(c), 1);
			}
			
			pList2.add(new Peak(100000, 0)); //we add a end-guard peak to end of the list
			                                 //makes subsquent iteration easier
			for(int c = 0; c < pList1.size(); c++){
				matchingGraph.addVertex(pList1.get(c), 2);
			}
			p1 = pList1.get(0);
			p2 = pList2.get(0);
			for(i = 0; i < pList1.size(); i++){
				p1 = pList1.get(i);
				//System.out.println("matching " + p1 + " starting " + j);
				for(; j < pList2.size(); j++){
					p2 = pList2.get(j);
					//scanning left to reach lower boundary for a theo-peak
					if(p1.getMass() - p2.getMass() > tolerance){
						continue;
					//once we surpass lower boundary, we check whether we are over 
					//upper boundary if it is we skip this theo-peak
					}else if(p2.getMass() - p1.getMass() > tolerance){
						break;
					//now we are within tolerance boundary	
					}else{
						//matchingGraph.addEdge(p2, p1);
						p = j;
						do{
							matchingGraph.addEdge(p2, p1);
							//System.out.println("matching: " + p2 + " ~ " + p1);
							p2 = pList2.get(++p);
						}while(p2.getMass() - p1.getMass() < tolerance);
						break;
//						for(int k = 1; k <= MAX_PROFILE_SPAN && j+k < pList2.size(); k++){
//							diff = pList2.get(j+k).getMass() - pList2.get(j+k-1).getMass() - 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
//							if(Math.abs(diff) < tolerance){
//								//System.out.println(actual.spectrumName + "\tmatched-iso-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j+k).getIntensity());
//								matchingGraph.addEdge(pList2.get(j+k), pList1.get(i));
//							}else{
//								break;
//							}
//						}
					}
					
				}
			}
			pList2.remove(pList2.size()-1);
			//System.out.println("there are total of " + matchingGraph.edgeSet().size() + " mathcing within tolerance");
			return matchingGraph;

	}
	
	public static List<Peak> getAnnotatedPeak(Spectrum s, TheoreticalSpectrum t, double tolerance){
		SimpleMatchingGraph g = t.getMatchGraph(s, tolerance);
		List<Peak> toBeRemoved = new ArrayList<Peak>();
		for(Iterator<Peak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak current = it.next();
			if(g.getNeighbors(current).size() > 0){
				toBeRemoved.add(current);
			}
		}
		Collections.sort(toBeRemoved, PeakIntensityComparator.comparator);
		return toBeRemoved;
	}
	
	public static void removeAnnotatedPeaks(Spectrum s, TheoreticalSpectrum t, double tolerance){
		removeAnnotatedPeaks(s, t, tolerance, 1000000);
	}
	
	public static void removeAnnotatedPeaks(Spectrum s, TheoreticalSpectrum t, double tolerance, int top){
		List<Peak> toBeRemoved = getAnnotatedPeak(s, t, tolerance);
		double matchedInt  = 0.0;
		double totalInt = 0.0;
		for(int i = 0; i < toBeRemoved.size(); i++){
			matchedInt += toBeRemoved.get(i).getIntensity();
		}
		for(int i = 0; i < s.getPeak().size(); i++){
			totalInt += s.getPeak().get(i).getIntensity();
		}
		System.out.println("scan: " + s.scanNumber + "\t" + t.peptide + "\texplainedInt: " + (matchedInt/totalInt));
		System.out.println("number of peaks removed " + toBeRemoved);
		System.out.println("before: " + s.getPeak().size());
		int begin = toBeRemoved.size() - top;
		begin = begin > 0 ? begin : 0;
		s.getPeaks().removeAll(toBeRemoved.subList(begin, toBeRemoved.size()));
		System.out.println("after: " + s.getPeak().size());
	}
	
	public static List<Peak> getTopPeaks(int N, List<Peak> pList){
		Collections.sort(pList, PeakIntensityComparator.comparator);
		int begin = pList.size() - N;
		begin = begin < 0 ? 0 : begin;
		return pList.subList(begin, pList.size());
	}
	
	public static void getUniqueSpectrum(String file){
		SpectrumLib lib = new SpectrumLib(file, "MGF");
		lib.printStat(SpectrumLib.NODETAIL);
		SpectrumLib lib2 = lib.Divide();
		lib2.printLibToFile("test.mgf", lib2);
	}
	
	
	public static void removeSUMOPeaks(Spectrum s){
		Peptide peptide = new Peptide("QQQTGG", 1);
		List<Peptide> pepList = new ArrayList<Peptide>();
		pepList.add(peptide);
		List<Peptide> linkedPeps = LinkedPeakScoreLearner.generateLinkedPeptides(pepList, s);
		Peptide linkPep = linkedPeps.get(1);
		TheoreticalSpectrum t = new TheoreticalSpectrum(linkPep, s.charge);
		SpectrumUtil.removeAnnotatedPeaks(s, t, 0.3);
	}
	
	
	public static void getSUMOPeaks(Spectrum s){
		MZXMLReader reader  = new MZXMLReader("../mixture_linked/linked_peptide_library/sumo_lib/20101008_Sumo_Library_4349_Bo.mzXML");
		s = reader.getSpectrum(4611);
		Peptide peptide = new Peptide("QQQTGG", 1);
		peptide.insertPTM(1, -17);
		List<Peptide> pepList = new ArrayList<Peptide>();
		pepList.add(peptide);
		List<Peptide> linkedPeps = LinkedPeakScoreLearner.generateLinkedPeptides(pepList, s, 'G');
		Peptide linkPep = linkedPeps.get(1);
		TheoreticalSpectrum t = new TheoreticalSpectrum(linkPep, s.charge);
		System.out.println("pep: " + linkPep);
		System.out.println("spec: " + s.parentMass + "\t" + s.charge);
		List<Peak> pList = SpectrumUtil.getAnnotatedPeak(s, t, 0.5);
		Collections.sort(pList, PeakMassComparator.comparator);
		s.setPeaks(pList);
		List<Spectrum> specList = new ArrayList();
		specList.add(s);
		SpectrumUtil.printSpectraToFile("../mixture_linked/SUMO.mgf", specList);
	}
	
	public static void getSUMORemovedSpectrum(){
		String spectrumFile = "..\\mixture_linked/linked_peptide_library/sumo_lib/human_sumo/Veronica_sumo_enrich/20120223/20120223_ananiav_SUMOscan_chymo_CID35_top15.mzXML";
		String annotationFile = "../mixture_linked/testAnnotation.txt";
		MZXMLReader reader  = new MZXMLReader(spectrumFile);
		List<String> lines = Utils.FileIOUtils.createListFromFile(annotationFile);
		List<Spectrum> specList = new ArrayList();
		for(int i = 0; i < lines.size(); i++){
			String[] tokens = lines.get(i).split("\\s+");
			Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[3]));
			String pep = tokens[13];
			System.out.println("tag is: " + pep);
			int ind = pep.lastIndexOf('+');
			pep = pep.substring(0, ind);
			System.out.println("tag is: " + pep);
			Peptide peptide = new Peptide(pep, 1);
			List<Peptide> pepList = new ArrayList<Peptide>();
			pepList.add(peptide);
			List<Peptide> linkedPeps = LinkedPeakScoreLearner.generateLinkedPeptides(pepList, s, 'G');
			Peptide linkPep = linkedPeps.get(1);
			System.out.println("pep: " + linkPep);
			System.out.println("spec: " + s.parentMass + "\t" + s.charge);
			TheoreticalSpectrum t = new TheoreticalSpectrum(linkPep, s.charge);
			SpectrumUtil.removeAnnotatedPeaks(s, t, 0.3);
			specList.add(s);
		}
		SpectrumUtil.printSpectraToFile("../mixture_linked/SUMORemoved.mgf", specList);

	}
	
	public static void generateMixSpectra(){
		String libfile1 = "..\\MSPLib\\Lib\\ecoli.msp";
		String libfile2 = "..\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib1 = new SpectrumLib(libfile1, "MSP");
		SpectrumLib lib2 = new SpectrumLib(libfile2, "MSP");
		//String libfile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		//SpectrumLib lib2 = new SpectrumLib(libfile2, "MGF");
		lib1.removeModSpectra();
		lib2.removeModSpectra();
		//lib2.createMix("..\\mixture_linked\\mixtures.mgf", 10000, 0.1, 0.00001, 1.0, 3, false, 5);
		for(int i = 0; i < 30; i++){
			SpectrumLib mixture = lib1.createRandomMix(lib2, 1000, 0.1, 0.001, 0.6, 2.5, false);
			System.out.println("Generated library size: " + mixture.getAllSpectrums().size());
			mixture.printLibToFile("..\\mixture_linked\\mixture.mgf_part"+(i+1), mixture);
		}
	}
	
	//compute the fraction of intensity from unfragmented precursors
	public static void precursorFraction(){
		String spectrumFile = "../mixture_linked/msdata/linked_peptide_library/ACG_disulfide_library/Orbi_Elite/PepLib1_300ng_trp_Elite_CID_35rep.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		while(reader.hasNext()){
			Spectrum s = reader.next();
			List<Peak> pList = s.getPeak();
			double totalInt = 0.0;
			double precursorInt = 0.0;
			for(int i = 0; i < pList.size(); i++){
				Peak p = pList.get(i);
				if(Math.abs(p.getMass() - s.parentMass) < 0.5){
					//System.out.println("matched " + p);
					precursorInt += p.getIntensity();
				}
				totalInt += p.getIntensity();
			}
			System.out.println(s.spectrumName +  "\t" + s.parentMass + "\t" + s.charge +"\tunfragmented intensity:\t" 
					+ (precursorInt/totalInt) + "\t" + precursorInt + "\t" + totalInt +"\t" + s.score);
		}
	}
	
	public static int getSNR(Spectrum s, double tolerance){
		Peptide pep = new Peptide(s.peptide, s.charge);
		String[] prefixIons = {"b", "b-H20", "b-NH3"};//, "b(iso)"};//"a", "a-H20", "a-NH3"};
		String[] suffixIons = {"y", "y-H20", "y-NH3"};//, "y(iso)"};
		TheoreticalSpectrum.prefixIons = prefixIons;
		TheoreticalSpectrum.suffixIons = suffixIons;
		TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
		TheoreticalSpectrum.addIsotopicPeaks(t, 1);
		s = new Spectrum(s); //make a copy of the spectrum because there is some processing below that will alter the spectrum
		s.removePrecursors(0.1);
		SimpleMatchingGraph g = t.getMatchGraph(s, tolerance);
		List<Peak> annotated = new ArrayList();
		List<Peak> unAnnotated = new ArrayList();
		List<Peak> allPeaks = new ArrayList();
		allPeaks.addAll(s.getPeak());
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				unAnnotated.add(p);
			}else{
				annotated.add(p);
			}
		}
		Collections.sort(allPeaks, PeakIntensityComparator.comparator);
		Collections.sort(annotated, PeakIntensityComparator.comparator);
		Collections.sort(unAnnotated, PeakIntensityComparator.comparator);
		double mean = 0.0;
		double topFract = 0.01;
		double minZ = 220;
		int length = unAnnotated.size();
		int count=0;
		for(int i = length-1; i > (int)(length*topFract); i--){
			mean+=unAnnotated.get(i).getIntensity();
			count++;
		}
		
		mean = mean / count;
		mean = allPeaks.get((int)(allPeaks.size()*0.25)).getIntensity();
		s.mergePeaks(s, 0.05);
		g = t.getMatchGraph(s, tolerance);
		unAnnotated.clear();
		annotated.clear();
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				unAnnotated.add(p);
			}else{
				annotated.add(p);
			}
		}
		
		int bigSignal = 0;
		for(int i = 0; i < annotated.size(); i++){
			double z = (annotated.get(i).getIntensity() / mean)*100;
			//annotated.get(i).setRank((int)z);
			if(z > minZ){
				bigSignal++;
				//System.out.println("big signal: " + annotated.get(i));
			}
		}
		
		int bigNoise = 0;
		for(int i = 0; i < unAnnotated.size(); i++){
			double z = (unAnnotated.get(i).getIntensity() / mean)*100;
			//unAnnotated.get(i).setRank((int)z);
			if(z > minZ){
				bigNoise++;
				//System.out.println("big nosie: " + unAnnotated.get(i));
			}
		}

		double medianNoi = allPeaks.get((int)(allPeaks.size()*0.25)).getIntensity();
		//double medianNoi = unAnnotated.get((int)(unAnnotated.size()*0.5)).getIntensity();
		if(annotated.size() == 0){
			System.out.println("warining: ");
			return 0;
		}
		double medianSig = annotated.get((int)(annotated.size()/2.0)).getIntensity();
		double totalInt = s.sumMagnitude();
		System.out.println(s.spectrumName + "\t" + s.peptide + "\t" + s.charge + "\ttotal peaks:\t" 
				+ annotated.size() + "\t" + unAnnotated.size() + "\tmean-noise level:\t" 
				+ mean + "\tmedian:\t" + medianNoi + "\tabove-" +minZ + "%:\t" + bigSignal + "\t" + bigNoise + "\ts2n-ratio:\t" + (medianSig/medianNoi)
				+ "\tabundance:\t" + s.score +"\t" + s.upperBound + "\t" + totalInt);
		return bigSignal;
	}
	
	public static void computeS2NR(){
		String spectrumDir = "../mixture_linked/msdata/EIF";
		String annotationFile = "..\\mixture_linked/test.txt";
		String outfile = "..\\mixture_linked\\test.mgf";
		List<Spectrum> specList = new ArrayList<Spectrum>();
		Map<String, MZXMLReader> readers = new HashMap<String, MZXMLReader>();
		try{
			BufferedReader buff = new BufferedReader(new FileReader(annotationFile));
			BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
			String line;
			String prevFile= "";
			MZXMLReader reader = null;
			int counter = 1; //one base index
			while((line=buff.readLine()) != null){
				String[] tokens = line.split("\\t");
				String file = tokens[0];
				file = file.replaceAll("\\s+", "");
				if(!readers.containsKey(file)){
					MZXMLReader newreader = new MZXMLReader(spectrumDir+"\\" + file);
					readers.put(file, newreader);
				}
				reader = readers.get(file);
				System.out.println("line\t"+line);
				Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[1]));
				String peptide = tokens[7];
				peptide = peptide.substring(2, peptide.length()-2);
				Peptide p = new Peptide(peptide, Integer.parseInt(tokens[6]));
				s.peptide=peptide;
				s.charge = Integer.parseInt(tokens[6]);
				s.scanNumber = counter++;
				int signals = SpectrumUtil.getSNR(s, 0.05);
				System.out.println(line + "\tsignals\t" + signals);
				//s.spectrumName = s.spectrumName + " PROTEIN: " + tokens[8]; 
				//out.write(s.toString());
				//out.write("\n");
				if(counter > 50){
					//break;
				}
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		//printSpectraToFile(outfile, specList);
	}
	
	
	public static void main(String[] args){
		//printLibPeptides();
		//analyzeMixtureSpectralMatches();
		//String filename = ".\\MSPLib\\Lib\\yeast.msp";
		///convertMSPToMGF(filename);
		//annotateSpectrumLib();
		//annotateSpectrumLibFromMzXML();
		if(Integer.parseInt(args[3]) == 7){
			annotateSpectrumLibFromMzXMLs(args[0], args[1], args[2]);
		}
		//annotateSpectrumLibFromMGF();
		//removeAnnotateFromMzXMLs();
		//getUniqueSpectrum("../mixture_linked/MSGFDBSearches/human_heck_e10_annotated_spectra.mgf");
		//insertParentMassInfo();
		//getSUMOPeaks(null);
		//getSUMORemovedSpectrum();
		//generateMixSpectra();
		//precursorFraction();
		//computeS2NR();
	}
}
