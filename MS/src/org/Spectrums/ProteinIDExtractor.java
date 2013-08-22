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
	String seqStr;
	Set<String> proteinIDs;
	Set<String> peptideIDs;
	Map<String, List<String>> peptideMap;
	Map<String, List<String>> proteinMap;
	Map<String, List<Object>> positionMap;
	
	private int pepIDIndex = 7;
	private int pepIDIndex2 = 8;
	
	public ProteinIDExtractor(String dbFile, String result){
		this.proteinDBFile = dbFile;
		this.searchResultFile = result;
		this.seqDB = new FastaSequence(this.proteinDBFile);
		this.seqStr = seqDB.getSubsequence(0, seqDB.getSize());
		System.out.println("resultFile: " + this.searchResultFile);
		this.outputFile = searchResultFile.split("\\.txt")[0]+".fasta";
		System.out.println("out: " + this.outputFile);
		parseResultFile();
		getPeptideProteinMap();
		//printProteins();
	}
	
	private void parseResultFile(){
		List<String> results = Utils.FileIOUtils.createListFromFile(this.searchResultFile);
		this.peptideIDs = new HashSet();
		for(int i = 0; i < results.size(); i++){
			String line = results.get(i);
			//System.out.println("line is: " + line);
			String[] tokens = line.split("\\t");
			String peptide = tokens[pepIDIndex];
			peptide = peptide.replaceAll("[0-9\\+\\-\\.\\_]", "");
			//peptide = peptide.substring(1, peptide.length()-1);
			//System.out.println("IDs is  : " + peptide);
			peptideIDs.add(peptide);
			String peptide2 = tokens[pepIDIndex2];
			peptide2 = peptide2.replaceAll("[0-9\\+\\-\\.\\_]", "");
			//peptideIDs.add(peptide2);
		}
	}
		
	private void getPeptideProteinMap(){
		this.peptideMap = new HashMap<String, List<String>>();
		this.proteinMap = new HashMap<String, List<String>>();
		this.positionMap = new HashMap<String, List<Object>>();
		this.proteinIDs = new HashSet();
		for(Iterator it = peptideIDs.iterator(); it.hasNext();){
			String id = (String)it.next();
			List<String> proteins = new ArrayList<String>();
			List<Object> positions = new ArrayList<Object>();
			this.peptideMap.put(id, proteins);
			this.positionMap.put(id, positions);
			//System.out.println("peptide IDs is: " + id);
			int index = this.seqStr.indexOf(id);
			long start = this.seqDB.getStartPosition(index);
			//System.out.println("starting residue: " + seqStr.charAt((int) start));
			while(index > 0){
				String protein = this.seqDB.getAnnotation(index);
				proteins.add(protein);
				positions.add(protein);
				positions.add((int)(index-start-1));
				List<String> peptides = null;
				if(this.proteinMap.containsKey(protein)){
					peptides = this.proteinMap.get(protein);
				}else{
					peptides = new ArrayList<String>();
					this.proteinMap.put(protein, peptides);
				}
				peptides.add(protein);
				//System.out.println("Protein ID is: "  + protein);
				proteinIDs.add(protein);
				index = this.seqStr.indexOf(id, index+1);
			}
		}
		System.out.println("After constructing positionmap size: " + this.positionMap.keySet().size());
		this.proteinIDs = this.proteinMap.keySet();
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
		List<String> results = Utils.FileIOUtils.createListFromFile(this.searchResultFile);
		String proteinIDPattern = "";
		int matchCount = 0;
		for(int i = 0; i < results.size(); i++){
			String line = results.get(i);
			String[] tokens = line.split("\\t");
			String peptide = tokens[pepIDIndex];
			peptide = peptide.replaceAll("[0-9\\+\\-\\.\\_]", "");
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
		ProteinIDExtractor IDS = new ProteinIDExtractor("../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta"
				, "../mixture_linked/test.txt");
		IDS.getPeptideReport();
		IDS.getProteinReport();
		IDS.printProteins(false);
		//IDS.getPeptidePairReport();
		//extractProteinsFromResultFiles("../mixture_linked/");
	}
}
