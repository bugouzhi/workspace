package org.Spectrums;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import IO.MZXMLReader;

/**
 * A class represent search result of a MS spectrum file
 * @author Jian
 *
 */
public class SearchResult{
	String resultFile;
	String spectrumFile;
	String databaseFile;
	String libraryfile;
	String separator = "\t";
	MZXMLReader reader ;
	SpectrumLib lib;
	public Map<Integer, AnnotatedSpectrum> resultMap;
	public Map<String, AnnotatedSpectrum> peptideMap;
	public Map<String, AnnotatedSpectrum> proteinMap;
	
	public SearchResult(String resultFile){
		this.resultFile = resultFile;
	}
	
	public SearchResult(String resultFile, String databaseFile, String libraryFile, String spectrumFile){
		this.resultFile = resultFile;
		this.databaseFile = databaseFile;
		this.libraryfile = libraryFile;
		this.spectrumFile = spectrumFile;
		this.lib = new SpectrumLib(this.libraryfile, "MGF");
		this.reader = new MZXMLReader(this.spectrumFile);
	}
	
	public SearchResult(String resultFile, String databaseFile, SpectrumLib lib, MZXMLReader reader){
		this.resultFile = resultFile;
		this.databaseFile = databaseFile;
		this.spectrumFile = spectrumFile;
		this.lib = lib;
		this.reader = reader;
	}
	
	public void parseResult(int fileInd, int scanInd, int pepInd, int chargeInd, int protInd, int scoreInd){
		Collection<String> results = Utils.FileIOUtils.createListFromFile(this.resultFile);
		this.resultMap = new HashMap();
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String result = it.next();
			String[] tokens = result.split(this.separator);
			//System.out.println("line is: " + result);
			if(tokens.length < scoreInd){
				continue;
			}
			int scan = Integer.parseInt(tokens[scanInd]);
			String pep = tokens[pepInd]+"." + tokens[chargeInd];
			String prot = tokens[protInd];
			double score = Double.parseDouble(tokens[scoreInd]);
			if(resultMap.containsKey(scan)){
				AnnotatedSpectrum s = resultMap.get(scan);	
				List<String> matches = (List<String>)s.getAnnotation().get("Matches");
				//System.out.println(matches);
				matches.add(pep);
				resultMap.put(s.scanNumber, s);
			}else{
				AnnotatedSpectrum s = new AnnotatedSpectrum();
				s.scanNumber = scan;
				s.peptide = pep;
				s.protein = prot;
				s.score = score;
				List<String> matches = new ArrayList<String>();
				matches.add(s.peptide);
				s.getAnnotation().put("Matches", matches);
				resultMap.put(s.scanNumber, s);
			}
		}
	}
	
	public Iterator<AnnotatedSpectrum> getResultIterator(){
		return this.resultMap.values().iterator();
	}
	
	public MixSSM getMatches(int scan){
		if(!this.resultMap.containsKey(scan)){
			System.out.println("warning cannot find result for spectrumScan " + scan);
			return null;
		}
		return MixSSM.createMixSSM(this.reader.getSpectrum(scan), 
				this.lib, 
				(List<String>)this.resultMap.get(scan).getAnnotation().get("Matches"));
	}
	
	
}
