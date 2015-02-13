package org.Spectrums;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import Utils.RKLookup;

import sequences.*;
/**
 * Given a sequence database and a list of search results from
 * database search engine:
 * - Extract list of protein IDs from a list of PSMs
 * - Mapped the peptide IDs to protein positions
 * * @author Jian Wang
 *
 */
public class ProteinIDExtractor {
	private String searchResultFile;
	private String proteinDBFile;
	private String outputFile;
	private BufferedWriter out;
	FastaSequence seqDB;
	char[] seqArray;
	String seqStr;
	Set<String> proteinIDs;
	Set<String> peptideIDs;
	RKLookup strMatches;
	Map<String, List<String>> peptideMap;
	Map<String, List<String>> proteinMap;
	Map<String, List<Object>> positionMap;
	
	private int pepIDIndex = 1;
	private int pepIDIndex2 = 1;
	
	
	
	public Map<String, List<String>> getPeptideMap() {
		return peptideMap;
	}

	public Map<String, List<String>> getProteinMap() {
		return proteinMap;
	}


	public Map<String, List<Object>> getPositionMap() {
		return positionMap;
	}


	public ProteinIDExtractor(String dbFile, String result){
		this.proteinDBFile = dbFile;
		this.searchResultFile = result;
		System.out.println("resultFile: " + this.searchResultFile);
		this.outputFile = searchResultFile.split("\\.txt")[0]+".fasta";
		System.out.println("out: " + this.outputFile);
		init();
		parseResultFile();
		getPeptideProteinMap();
		//printProteins();
	}
	
	public ProteinIDExtractor(String dbFile, String result, int pepInd){
		this.proteinDBFile = dbFile;
		this.searchResultFile = result;
		this.pepIDIndex = pepInd;
		System.out.println("resultFile: " + this.searchResultFile);
		this.outputFile = searchResultFile.split("\\.txt")[0]+".fasta";
		System.out.println("out: " + this.outputFile);
		init();
		parseResultFile();
		getPeptideProteinMap();
		//printProteins();
	}
	
	public ProteinIDExtractor(List IDs, String dbFile){
		this.peptideIDs = new HashSet();
		//System.out.println("Protein fasta: " + dbFile);
		this.proteinDBFile = dbFile;
		this.seqDB = new FastaSequence(this.proteinDBFile);
		this.seqStr = seqDB.getSubsequence(0, seqDB.getSize());
		System.out.println("resultFile: " + this.searchResultFile);
		//this.outputFile = searchResultFile.split("\\.txt")[0]+".fasta";
		//System.out.println("out: " + this.outputFile);
		//System.out.println("IDs-list size: " + IDs.size());
		for(int i = 0; i < IDs.size(); i++){
			Spectrum s = (Spectrum)IDs.get(i);
			String peptide = s.peptide;
			//System.out.println("peptide: " + peptide);
			//peptide = getStrippedSeq(peptide);
			peptideIDs.add(peptide);
		}
		init();
		getPeptideProteinMap();
		//getPeptideReport();
	}
	
	
	public ProteinIDExtractor(Set<String> IDs, String dbFile){
		this.peptideIDs = new HashSet();
		//System.out.println("Protein fasta: " + dbFile);
		this.proteinDBFile = dbFile;
		this.seqDB = new FastaSequence(this.proteinDBFile);
		this.seqStr = seqDB.getSubsequence(0, seqDB.getSize());
		System.out.println("resultFile: " + this.searchResultFile);
		//this.outputFile = searchResultFile.split("\\.txt")[0]+".fasta";
		//System.out.println("out: " + this.outputFile);
		//System.out.println("IDs-list size: " + IDs.size());
		for(Iterator<String> it = IDs.iterator(); it.hasNext();){
			this.peptideIDs.add(it.next());
		}
		init();
		getPeptideProteinMap();
		//getPeptideReport();
	}
	
	
	private void init(){
		this.seqDB = new FastaSequence(this.proteinDBFile);
		this.seqStr = seqDB.getSubsequence(0, seqDB.getSize());
		this.seqArray = seqStr.toCharArray();
	}
	
	private String getStrippedSeq(String pep){
		return Utils.StringUtils.getStrippedSeq(pep);
	}
	
	private void parseResultFile(){
		List<String> results = Utils.FileIOUtils.createListFromFile(this.searchResultFile);
		this.peptideIDs = new HashSet();
		for(int i = 0; i < results.size(); i++){
			String line = results.get(i);
			//System.out.println("line is: " + line);
			String[] tokens = line.split("\\t");
			if(tokens.length < this.pepIDIndex || line.startsWith("#") || line.startsWith("SpectrumFile")){
				continue;
			}
			String peptide = tokens[pepIDIndex];
			if(peptide.contains("!")){
				String[] peptides = peptide.split("!");
				for(int p = 0; p < peptides.length; p++){
					peptide = Utils.StringUtils.getStrippedSeq(peptides[p]);
					if(peptides.length > 1){
						//System.out.println("adding peptide: " + peptide);
					}
					peptideIDs.add(peptide);
				}
			}else{
				peptide = Utils.StringUtils.getStrippedSeq(peptide);
				//peptide = peptide.substring(1, peptide.length()-1);
				//System.out.println("IDs is  : " + peptide);
				peptideIDs.add(peptide);
				//String peptide2 = tokens[pepIDIndex2];
				//peptide2 = peptide2.replaceAll("[0-9\\+\\-\\.\\_\\(\\)]", "");
				//peptideIDs.add(peptide2);
			}
		}
		//System.out.println("Done parsing results");
	}
		
	private void getPeptideProteinMap(){
		this.peptideMap = new HashMap<String, List<String>>();
		this.proteinMap = new HashMap<String, List<String>>();
		this.positionMap = new HashMap<String, List<Object>>();
		this.proteinIDs = new HashSet();
		int counter =0;
		Map<String, List<String>> modSeqMap = new HashMap<String, List<String>>();
		for(Iterator<String> it = this.peptideIDs.iterator(); it.hasNext();){
			String pep = it.next();
			String stripped = getStrippedSeq(pep);
			if(modSeqMap.containsKey(stripped)){
				modSeqMap.get(stripped).add(pep);
			}else{
				List<String> mods = new ArrayList();
				mods.add(pep);
				modSeqMap.put(stripped, mods);
			}
		}
		RKLookup strMatcher = new RKLookup(modSeqMap.keySet(), 5);
		Map<String, List<Integer>> matches = strMatcher.matches(this.seqStr);
		for(Iterator<String> it = matches.keySet().iterator(); it.hasNext();){
			String pep = it.next();
			List<Integer> positions = matches.get(pep);
			List<String> proteins = new ArrayList<String>();
			List<Object> pos = new ArrayList<Object>();
			//System.out.println("Peptide: " + pep +"\tmatches: " + positions.size());
			for(int i = 0; i < positions.size(); i++){
				int ind = positions.get(i);
				String prot = this.seqDB.getAnnotation(ind);
				List<String> peps;
				if(this.proteinMap.containsKey(prot)){
					peps = this.proteinMap.get(prot);
				}else{
					peps = new ArrayList<String>();
					this.proteinMap.put(prot, peps);
				}
				proteins.add(prot);
				pos.add(prot);
				pos.add(ind-this.seqDB.getStartPosition(ind-1));
			}
			List<String> modSeq = modSeqMap.get(pep);
			for(int i = 0; i < modSeq.size(); i++){
				for(int j = 0; j < proteins.size(); j++){
					this.proteinMap.get(proteins.get(j)).add(modSeq.get(i));
				}
				//System.out.println("putting mod " + modSeq.get(i));
				this.peptideMap.put(modSeq.get(i), proteins);
				this.positionMap.put(modSeq.get(i), pos);
			}
		}
		System.out.println("checking pep2\t" + this.peptideMap.containsKey("+42.011MQNDAGEFVDLYVPR"));
		
		System.out.println("After constructing positionmap size: " + this.positionMap.keySet().size());
		this.proteinIDs = this.proteinMap.keySet();
	}
	
	
	//get a list of peptides not shared by proteins
	public Set<String> getNonSharedPeps(){
		Set<String> nonDegenerate = new HashSet<String>();
		System.out.println("size: " + this.peptideMap.keySet().size() + "\t" + this.peptideIDs.size());
		for(Iterator<String> it = this.peptideMap.keySet().iterator(); it.hasNext();){
			String pep = it.next();
			List<String> protIds = this.peptideMap.get(pep);
			if(protIds.size() == 1){
				nonDegenerate.add(pep);
			}
		}
		System.out.println("Total non-degenerate " + nonDegenerate.size());
		return nonDegenerate;
	}
	
	public Set<String> getNonSharedPeps(String prot){
		Set<String> nonDegenerate = new HashSet<String>();
		List<String> peps = this.proteinMap.get(prot);
		if(peps == null){
			return nonDegenerate;
		}
		for(Iterator<String> it = peps.iterator(); it.hasNext();){
			String pep = it.next();
			List<String> protIds = this.peptideMap.get(pep);
			if(protIds.size() == 1){
				nonDegenerate.add(pep);
				//System.out.println("Protein " + prot +"\tuniq-pep:\t" + pep);
			}
		}
		//System.out.println("Total non-degenerate " + nonDegenerate.size());
		return nonDegenerate;
	}
	
	//get a list of proteins do not contain shared peptides
	public Set<String> getProtWithNoSharePeps(){
		Set<String> prots =  new HashSet<String>();
		Set<String> peps = getNonSharedPeps();
		for(Iterator<String> it = peps.iterator(); it.hasNext();){
			String pep = it.next();
			prots.add(this.peptideMap.get(pep).get(0));
		}
		return prots;
	}
	
	
	public void printProteins(boolean printDecoy){
		try{
			this.out = new BufferedWriter(new FileWriter(this.outputFile));
			for(Iterator it = proteinIDs.iterator(); it.hasNext();){
				String protein = (String)it.next();
				String protStr = this.seqDB.getMatchingEntry(protein);
				if(printDecoy || !protein.startsWith("X_") ){
					out.append(">"+protein+"\n");
					out.append(protStr +"\n");
				}
			}
			//print decoy
			for(Iterator it = proteinIDs.iterator(); it.hasNext();){
				String protein = (String)it.next();
				String protStr = this.seqDB.getMatchingEntry(protein);
				if(!protein.startsWith("X_") && printDecoy){
					out.append(">X_"+protein +"\n");
					out.append(new StringBuffer(protStr).reverse().toString() +"\n");
				}
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public void getPeptideReport(){
		String proteinIDPattern = "";
		int matchCount = 0;
		for(Iterator<String> it = this.peptideMap.keySet().iterator(); it.hasNext();){
			String peptide = it.next();
			//peptide = peptide.substring(1, peptide.length()-1);
			List<String> proteins = this.peptideMap.get(peptide);
			System.out.print("ID: " + peptide + "\t");
			for(int j = 0; j < proteins.size(); j++){
				System.out.print(proteins.get(j) + ";");
			}
			System.out.println();
		}
	}
	
	public void getProteinReport(){
		for(Iterator<String> it = this.proteinMap.keySet().iterator(); it.hasNext();){
			String protID = it.next();
			List<String> pepList = this.proteinMap.get(protID);
			System.out.println("Protein: " + protID + "\tmapped peps:\t" + pepList.size());
		}
	}
	
	public void getPeptidePairReport(){
		List<String> results = Utils.FileIOUtils.createListFromFile(this.searchResultFile);
		for(int i = 0; i < results.size(); i++){
			String line = results.get(i);
			String[] tokens = line.split("\\s+");
			String peptide1 = tokens[this.pepIDIndex];
			String peptide2 = tokens[this.pepIDIndex2];
			peptide1 = peptide1.replaceAll("[0-9\\+\\-\\.r]", "");
			peptide2 = peptide2.replaceAll("[0-9\\+\\-\\.r]", "");
			//peptide = peptide.substring(1, peptide.length()-1);
			List<String> proteins1 = this.peptideMap.get(peptide1);
			System.out.print("ID1: " + peptide1 + "\t");
			for(int j = 0; j < proteins1.size(); j++){
				String prot = proteins1.get(j);
				//prot = prot.split("\\s+")[0];
				System.out.print(prot + ";");
			}
			List<String> proteins2 = this.peptideMap.get(peptide2);
			System.out.print("\t");
			System.out.print("ID2: " + peptide2 + "\t");
			for(int j = 0; j < proteins2.size(); j++){
				String prot = proteins2.get(j);
				//prot = prot.split("\\s+")[0];
				System.out.print(prot + ";");
			}
			System.out.println();
		}
	}
	
	public Map<String, List<Object>> getPositonMap(){
		return this.positionMap;
	}
	
	//mapp ID to proteins
	public static void mapID(String ID, Map<String, List<Integer>> idMap, List<String> seqs){
		if(ID.startsWith("r")){
			ID = ID.substring(1);
		}
		if(idMap.containsKey(ID)){
			return;
		}
		idMap.put(ID, new ArrayList<Integer>());
		List<Integer> positions = idMap.get(ID);
		for(int j = 0; j < seqs.size(); j++){
			String chainSeq = seqs.get(j);
			//position of peptide in chain
			int index = chainSeq.indexOf(ID);
			//System.out.println("chain: " + j + " index: " + index);
			//linked position, now assume to be first lysine
			int position = index + ID.indexOf('K');
			if(index >=0){
				positions.add(j);
				positions.add(index);
				positions.add(position);
			}
		}
	}
	
	public static void extractProteinsFromResultFiles(String directory){
		File dir = new File(directory);
		File[] files = dir.listFiles();
		boolean printDecoy = false;
		for(int i =0; i < files.length; i++){
			File current = files[i];
			//System.out.println("file is: " + current.getName());
			String name = current.getName();
			if(name.matches("Rabbit_proteasome_uniprot_rabbitseqs_M13[0-9]_msgfdbIDs\\.txt")){
				System.out.println("file is: " + current.getName());
				ProteinIDExtractor IDS = new ProteinIDExtractor("../mixture_linked/database/rabbit_uniprot_allproteins_plusDecoy.fasta"
					, directory+"/"+current.getName());
				IDS.printProteins(printDecoy);
				
			}
		}
	}
	
	public static void main(String[] args){
		ProteinIDExtractor IDS = new ProteinIDExtractor("../mixture_linked/database/UPS2.fasta"
				, "../mixture_linked/18476_REP3_400fmol_UPS2_500ugHumanLysate_peakview - Peptides.txt");
		IDS.getPeptideReport();
		//IDS.getProteinReport();
		//IDS.printProteins(false);
		//IDS.getPeptidePairReport();
		//extractProteinsFromResultFiles("../mixture_linked/");
	}
}
