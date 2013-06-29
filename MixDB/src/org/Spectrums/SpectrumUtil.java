package org.Spectrums;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
				//LinkedPeptide lp = new LinkedPeptide(peps[0]+"--"+peps[1], s.charge,  3, 6);
				//s.setPeptide(p.getPeptide());//+"."+s.charge);
				//s.parentMass = p.getParentmass();
				//s.charge = p.getCharge();
				//s.parentMass = lp.getParentmass();
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
	
	public static void annotateSpectrumLib(){
		String libFile = "..\\mixture_linked\\yeast_data/klc_122007p_yeast_digest1.mgf";
		String annotationFile = "..\\mixture_linked\\yeast1_MSPLIT_result.txt";
		SpectrumLib lib1 = new SpectrumLib(libFile, "MGF");
		loadAnnotationFromFile(annotationFile, lib1);
		lib1.printLibToFile("..\\mixture_linked\\yeast1_annotatedLib.mgf", lib1);
	}
	
	public static void annotateSpectrumLibFromMzXML(){
		String libFile = "..\\mixture_linked\\human_heck_data\\data\\090121_NM_Trypsin_30.mzXML";
		String annotationFile = "..\\mixture_linked\\testAnnotation.txt";
		String outfile = "..\\mixture_linked\\test.mgf";
		MZXMLReader reader = new MZXMLReader(libFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		loadAnnotationFromFile(annotationFile, specList);
		printSpectraToFile(outfile, getAnnotatedSpectra(specList));
	}
	
	public static void annotateSpectrumLibFromMGF(){
		String libFile = "..\\mixture_linked\\MSGFDBSearches\\cid_090121_NM_Trypsin_30.mgf";
		String annotationFile = "..\\mixture_linked\\MSGFDBSearches\\annotation_30_probe10.txt";
		String outfile = "..\\mixture_linked\\test.mgf";
		LargeSpectrumLibIterator<Spectrum> it = new LargeSpectrumLibIterator(libFile);
		List<Spectrum> specList = it.readAllSpectra();
		loadAnnotationFromFile(annotationFile, specList);
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
		List<Peak> toBeRemoved = getAnnotatedPeak(s, t, tolerance);
		System.out.println("number of peaks removed " + toBeRemoved);
		System.out.println("before: " + s.getPeak().size());
		s.getPeaks().removeAll(toBeRemoved);
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
	
	public static void main(String[] args){
		//printLibPeptides();
		//analyzeMixtureSpectralMatches();
		//String filename = ".\\MSPLib\\Lib\\yeast.msp";
		///convertMSPToMGF(filename);
		//annotateSpectrumLib();
		annotateSpectrumLibFromMzXML();
		//annotateSpectrumLibFromMGF();
		//getUniqueSpectrum("../mixture_linked/MSGFDBSearches/human_heck_e10_annotated_spectra.mgf");
		//insertParentMassInfo();
	}
}
