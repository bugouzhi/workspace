package org.Spectrums;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import IO.MZXMLReader;
import Utils.ArrayUtils;

/**
 * Report search results related stats
 * @author Jian
 *
 */
public class SearchResultReporter {
	String spectrumFile = "";
	String resultsFile = "";
	String databaseFile = "";
	int spectFileInd = 0;
	int pepInd = 7;
	int protInd = 8;
	MZXMLReader reader;
	TDAStat tdaStat;
	
	
	public String getSpectrumFile() {
		return spectrumFile;
	}

	public void setSpectrumFile(String spectrumFile) {
		this.spectrumFile = spectrumFile;
	}

	public SearchResultReporter(String spectrumFile, String resultFile, int pepInd, int protInd, int scoreInd, int mode){
		//System.out.println(spectrumFile);
		System.out.println(resultFile);
		this.spectrumFile = spectrumFile;
		this.resultsFile = resultFile;
		this.tdaStat = new TDAStat(resultFile, pepInd, protInd, scoreInd, mode);
		File f = new File(spectrumFile);
		if(f.exists()){
		//	this.reader = new MZXMLReader(spectrumFile);
		}else{
			System.err.println("warining: cannot find spectrum file, running with spectrum information");
		}
	}
	
	public SearchResultReporter(String resultFile){
		this.spectrumFile = getSpectrumFile(resultFile, this.spectFileInd);
		this.resultsFile = resultFile;
		this.tdaStat = new TDAStat(resultFile);
		//this.reader = new MZXMLReader(spectrumFile);
	}
	
	public static String getSpectrumFile(String resultFile, int spectFileInd){
		List<String> lines = Utils.FileIOUtils.readLines(resultFile, 5);
		for(int i = lines.size()-1; i > 0; i--){
			String[] tokens = lines.get(i).split("\\t");
			if(tokens[0].startsWith("#") || tokens.length < 5){
				continue;
			}
			String putativePath = tokens[spectFileInd];
			String[] tokens2 = putativePath.split("[/\\\\]");
			return tokens2[tokens2.length-1];
		}
		return "";
	}
	
	public void getStatSummary(){
		System.out.println("Data file: " + this.spectrumFile);
		if(this.reader != null){
			System.out.println("Total Scans: " + this.reader.getParser().getScanCount());
	//		int[] scanCounts = reader.getSpectrumStat();
	//		System.out.println("Total MS/MS scans: " + scanCounts[1]);
			System.out.println("Total MS/MS scans: " + 100);
		}
		tdaStat.getSummary();
	}
	
	public Collection<Spectrum> getSubResults(SpectrumFilter[] filters){
		Collection<AnnotatedSpectrum> results = tdaStat.getPeptideResults();
		//Collection<AnnotatedSpectrum> results = tdaStat.getResults();
		Collection<Spectrum> filtered = new ArrayList<Spectrum>();
		for(Iterator<AnnotatedSpectrum> it = results.iterator(); it.hasNext();){
			while(it.hasNext()){
				Spectrum s = it.next();
				boolean pass = true;
				for(int i = 0; i < filters.length; i++){
					SpectrumFilter filter = filters[i];
					if(!filter.accept(s)){
						pass = false;
						continue;
					}
					//System.out.println("passing fitler " + i);
				}
				if(pass){
					filtered.add(s);
				}
			}
		}
		return filtered;
	}
	
	public static int getOverlapPeptide(Collection<Spectrum> results1, Collection<Spectrum> results2){
		Collection list = new ArrayList();
		for(Iterator<Spectrum> it = results1.iterator(); it.hasNext();){
			list.add(it.next().peptide);
		}
		
		for(Iterator<Spectrum> it = results2.iterator(); it.hasNext();){
			list.add(it.next().peptide);
		}
		Map<Object, Integer> pepMap = getCount(list);
		Map<Object, Integer> countMap = getCount(pepMap.values());
		return countMap.get(new Integer(2));
	}
	
	public static Map<Object, Integer> getOverlapPeptide(Collection<Collection<Spectrum>> results){
		Collection pepList = new ArrayList();
		Map<Object, Integer> combCountMap = new HashMap<Object, Integer>();
		for(Iterator<Collection<Spectrum>> it = results.iterator(); it.hasNext();){
			for(Iterator<Spectrum> it2 = it.next().iterator(); it2.hasNext();){
				pepList.add(it2.next().peptide);
			}
		}
		Map<Object, Integer> pepMap = getCount(pepList);
		Map<Object, Integer> countMap = getCount(pepMap.values());
		return countMap;
	}
	
	
	public static Map<Object, Integer> getCount(Collection manyObjects){
		Map<Object, Integer> objMap = new HashMap<Object, Integer>();
		for(Iterator it = manyObjects.iterator(); it.hasNext();){
			Object current = it.next();
			if(objMap.containsKey(current)){
				Integer count = objMap.get(current);
				count = count+1;
				objMap.put(current, count);
			}else{
				objMap.put(current, 1);
			}
		}
		
		return objMap;
	}
	

	
	public Object[] getCustomReport1(){
		SpectrumFilter unModFilter = new SpectrumPepFilter("[^0-9]+");
		SpectrumFilter upsProtFilter = new SpectrumProteinFilter(".*HUMAN.*");
		//SpectrumFilter upsProtFilter = new SpectrumProteinFilter(".*SV=[0-8]_[0\\.]*5[0]*$");
		SpectrumFilter upsAbundFilter1 = new SpectrumProteinFilter(".*HUMAN.*_50000$");
		SpectrumFilter upsAbundFilter2 = new SpectrumProteinFilter(".*HUMAN.*_5000$");
		SpectrumFilter upsAbundFilter3 = new SpectrumProteinFilter(".*HUMAN.*_500$");
		SpectrumFilter upsAbundFilter4 = new SpectrumProteinFilter(".*HUMAN.*_50$");
		SpectrumFilter upsAbundFilter5 = new SpectrumProteinFilter(".*HUMAN.*_5$");
		SpectrumFilter upsAbundFilter6 = new SpectrumProteinFilter(".*HUMAN.*_0.5$");
		//this.tdaStat.getSpecCount(0.01, 0.012);
		this.tdaStat.mapProtein("../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta");
		SpectrumFilter fdrFilter = new SpectrumFDRFilter(0.010002, SpectrumFDRFilter.PEP);
		SpectrumFilter[] filters = new SpectrumFilter[]{fdrFilter, unModFilter};
		SpectrumFilter[] filters2 = new SpectrumFilter[]{fdrFilter, unModFilter, upsProtFilter};
		// int[] scanCounts = reader.getSpectrumStat();
		int[] scanCounts = new int[]{69825, 67829, 66000};
		
		Collection<Collection> subResults = new ArrayList<Collection>();
		subResults.add(this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter}));
		subResults.add(this.getSubResults(filters));
		subResults.add(this.getSubResults(filters2));
		
		
		//Peptide IDs stat
		subResults.add(this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter, upsProtFilter, upsAbundFilter1}));
		subResults.add(this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter, upsProtFilter, upsAbundFilter2}));
		subResults.add(this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter, upsProtFilter, upsAbundFilter3}));
		subResults.add(this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter, upsProtFilter, upsAbundFilter4}));
		subResults.add(this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter, upsProtFilter, upsAbundFilter5}));
		subResults.add(this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter, upsProtFilter, upsAbundFilter6}));
		List<Spectrum> filtered = new ArrayList();
		filtered.addAll(subResults.iterator().next());
		
		ProteinIDExtractor protID = new ProteinIDExtractor(filtered, "../mixture_linked/database/Human_allproteins_withShortname.fasta");
		Set<String> peps = protID.getNonSharedPeps();
		//peptide stats
		for(int i = 0; i < filtered.size(); i++){
			AnnotatedSpectrum s = (AnnotatedSpectrum)filtered.get(i);
			int unique = 0;
			String protein = "";
			if(peps.contains(s.peptide)){
				unique = 1;
			}
			if(protID.peptideMap.containsKey(s.peptide)){
				protein = protID.peptideMap.get(s.peptide).get(0);
			}
			System.out.println("PepList\t" + this.resultsFile + "\t" + Utils.StringUtils.getPepSeq(s.peptide) + "\t"
					+ s.score + "\t" + s.getAnnotation().get("pepfdr") +"\t" + unique + "\t" + protein + "\t" + s.protein + "\t" + s.getAnnotation().get("specCount"));
		}
		//Set<String> peps = protID.peptideIDs;
		Map<String, AnnotatedSpectrum> proteinMap = this.tdaStat.proteinMap;
		Set<String> toBeRemove = new HashSet();
		Map<String, AnnotatedSpectrum> newMap = new HashMap<String, AnnotatedSpectrum>();
		int decoys = 0;
		for(Iterator<String> it = proteinMap.keySet().iterator(); it.hasNext();){
			String prot =it.next();
			String pep = proteinMap.get(prot).peptide;
			pep = Utils.StringUtils.getStrippedSeq(pep);
			AnnotatedSpectrum s = proteinMap.get(prot);
			if(prot.startsWith("DECOY") || prot.startsWith("X_")){
				s.getAnnotation().put("pepcount", 1);
				s.getAnnotation().put("uniquePepcount", 1);
				newMap.put(s.protein, s);
			}else{
				//System.out.println("Protein is: " + s.protein + "\tpep\t" + pep);
				if(protID.peptideMap.containsKey(pep)){
					//if(prot.contains("DECOY") || prot.startsWith("X_")){
					//	decoys++;
					//}
					//if(s.protein.contains("RS21")){
					//	System.out.println("before prot: " + s.protein);
					//}
					if(peps.contains(pep)){
						s.protein = protID.peptideMap.get(pep).get(0);
					}
					//System.out.println("prot: " + s.protein);
					//System.out.println("protein name: " + pname);
					//System.out.println("pep-count is: " + protID.proteinMap.get(s.protein).size());
					Set<String> uniques = protID.getNonSharedPeps(s.protein);
					int specCount = 0;
					//System.out.println("Protein is: " + s.protein + "\t" + uniques.size());
					if(uniques.size() > 0){
						s.getAnnotation().put("pepcount", protID.proteinMap.get(s.protein).size());
						s.getAnnotation().put("uniquePepcount", uniques.size());
						newMap.put(s.protein, s);
					}
				}
			}
		}
		System.out.println("removing: " + toBeRemove.size());
		//System.out.println("decoys " + decoys);
		//for(Iterator<String> it = toBeRemove.iterator(); it.hasNext();){
		//	proteinMap.remove(it.next());
		//}
		this.tdaStat.proteinMap = newMap;
		this.tdaStat.getProteinFDR();
		
		
		//preferential assignment peptide to UPS proteins
		for(int i = 0; i < this.tdaStat.getResults().size(); i++){
			Spectrum s = this.tdaStat.getResults().get(i);
			List<String> prots = protID.peptideMap.get(s.peptide);
			if(prots != null && prots.size() > 0){
				//System.out.println("mapping: " + prots.get(0));
				s.protein = prots.get(0);
				
			}
		}
		
		int[] pepNum = new int[]{0,1,2,3,4,5,6,7,8,9,1000};
		int[] pepCounts = new int[pepNum.length];
		for(Iterator<String> it = newMap.keySet().iterator(); it.hasNext();){
			String prot =it.next();
			String pep = newMap.get(prot).peptide; 
			AnnotatedSpectrum s = newMap.get(prot);
			//System.out.println(s.protein);
			int pepcount = (Integer)s.getAnnotation().get("pepcount");
			int uniquecount = (Integer)s.getAnnotation().get("uniquePepcount");
			double protfdr = (Double)s.getAnnotation().get("Protfdr");
			String shortProtName;
			if(s.protein.contains("RepID")){
				shortProtName = s.protein.split("RepID=")[1];
			}else{
				shortProtName = s.protein.split("_")[0];
			}
			System.out.println("ProtList\t" + this.resultsFile +"\t" + shortProtName + "\t" 
					+ pepcount + "\t" + uniquecount + "\t" + protfdr
					+ "\t" + s.getAnnotation().get("specCount"));
		
			if(protfdr <= 0.01){
				int ind = ArrayUtils.getIntervalIndex(pepcount, pepNum);
				pepCounts[ind]++;
			}
		}
		
		System.out.print("SpectrumFile\tScans\tIDs\tUnique-Peptides\tUnmod-Peptides\tUps-Peptides\n");
		System.out.print(spectrumFile + "\t");
		if(this.reader != null) {
			System.out.print(this.reader.getParser().getScanCount() +"\t");
		}else{
			System.out.print("0\t");
		}
		//Total MS2 spetra
		System.out.print(scanCounts[1] +"\t");
		//total PSMs
		System.out.print(this.tdaStat.getPSMs()[0] +"\t");
		//peptide-leve IDs
		for(Iterator<Collection> it = subResults.iterator(); it.hasNext();){
			System.out.print(it.next().size() +"\t");
		}
		for(int i = 0; i < this.tdaStat.Prots.length; i++){
			System.out.print(this.tdaStat.Prots[i] +"\t");
		}
		System.out.println();
		
		//pepcount
		 System.out.println(Arrays.toString(pepCounts));
		return subResults.toArray();
	}
	
	public void getPeptideReport(){
		SpectrumFilter fdrFilter = new SpectrumFDRFilter(0.012, SpectrumFDRFilter.PEP);
		SpectrumFilter unModFilter = new SpectrumPepFilter("[^0-9]+");
		Collection<Spectrum> subResult = this.getSubResults(new SpectrumFilter[]{fdrFilter, unModFilter});
		System.out.println("total result: " + subResult.size());
		//for(Iterator<Spectrum> it = subResult.iterator(); it.hasNext();){
		//	System.out.println(it.next());
		//}
	}
	
	public static void testSearchResultReport(){
		String spectFilePath = "../mixture_linked/msdata/";
		String resultFilePath = "../mixture_linked/SWATH/PeakView_withRtCalc/";
		File spectrumFiles = new File(spectFilePath);
		File resultFiles = new File(resultFilePath);
		File[] files = resultFiles.listFiles();
		for(int i = 0; i < files.length; i++){
			if(files[i].isDirectory() || !files[i].getName().matches("18511.*Peptides\\.txt")){
				continue;
			}
			String searchResult = files[i].toString();
			String spectrumFile = spectFilePath + getSpectrumFile(searchResult, 0);
			//SearchResultReporter reporter = new SearchResultReporter(spectrumFile, searchResult, 7,8,11,1);
			//SearchResultReporter reporter = new SearchResultReporter(spectrumFile, searchResult,4,8,32,-1);
			SearchResultReporter reporter = new SearchResultReporter(spectrumFile, searchResult,1,0,5,1);
			//reporter.getStatSummary();
			reporter.getCustomReport1();
		}
	}
	
	public static void testSearchResultReport(String file){
		file = "../mixture_linked//out.txt";
		String searchResult = file;
		//SearchResultReporter reporter = new SearchResultReporter(spectrumFile, searchResult, 7,8,11,1);
		SearchResultReporter reporter = new SearchResultReporter("", searchResult,4,8,32,-1);
		reporter.getPeptideReport();
	}
	
	
	public static void testSearchResultReportPair(){
		String spectFilePath = "../mixture_linked/msdata/UPS12_Human/";
		String resultFilePath = "../mixture_linked/UPS_EcoliREP2_REP3_searches/IDA/UPSEcoli_REP3_newStock2_msgfdb_0.txt ";
		File spectrumFiles = new File(spectFilePath);
		File resultFiles = new File(resultFilePath);
		File[] files = resultFiles.listFiles();
		List<Object[]> results = new ArrayList<Object[]>();
		for(int i = 0; i < files.length; i++){
			if(files[i].isDirectory()){
				continue;
			}
			String searchResult = files[i].toString();
			String spectrumFile = spectFilePath + getSpectrumFile(searchResult, 0);
			SearchResultReporter reporter = new SearchResultReporter(spectrumFile, searchResult, 7,8,11,1);
			//SearchResultReporter reporter = new SearchResultReporter(spectrumFile, searchResult,4,8,31,-1);

			Object[] result = reporter.getCustomReport1();
			results.add(result);
		}
		
		int groupsize = 3;
		for(int i = 0; i < results.size(); i+=groupsize){
			for(int j = 0;  j < groupsize; j++){
				for(int l = j+1; l < groupsize; l++){
					System.out.print("Two-overlapCount:\t");
					for(int k = 0; k < results.get(0).length; k++){
						Collection<Collection<Spectrum>> resultList = new ArrayList();	
						resultList.add((Collection)results.get(i+j)[k]);
						resultList.add((Collection)results.get(i+l)[k]);
						Map<Object, Integer> countMap = getOverlapPeptide(resultList);
						System.out.print(countMap.get(2) +"\t");
					}
					System.out.println();
				}
			}
		}

		for(int i = 0; i < results.size(); i+=groupsize){
			System.out.print("ThreeoverlapCount:\t");
			for(int k = 0; k < results.get(0).length; k++){
				Collection<Collection<Spectrum>> resultList = new ArrayList();	
				for(int j = 0;  j < groupsize; j++){
					resultList.add((Collection)results.get(i+j)[k]);
				}
				Map<Object, Integer> countMap = getOverlapPeptide(resultList);
				System.out.print(countMap.get(3) +"\t");
			}
			System.out.println();
		}
	}
	
	
	public static void main(String[] args){
		testSearchResultReport();
		//testSearchResultReport(args[0]);
		//testSearchResultReportPair();
	}
	
	
	
}
