package org.Spectrums;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import Utils.FileIOUtils;
import org.biojava.bio.structure.*;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.*;

/**
 * We map IDs to a structure to use strucural information
 * to validate our IDs (particularly crosslinked IDs)
 * @author Jian Wang
 *
 */
public class IDs2StructureMapper {
	private String resultFile;
	private String pdbID;
	private String pdbFile;
	private List<String> IDs1;
	private List<String> IDs2;
	private Structure structure;
	private int idIndex1 = 10;
	private int idIndex2 = 12;
	private List<String> seqs;
	private Map<String, List> idMap;

	
	public IDs2StructureMapper(String resultFile, String pdbFile){
		this.resultFile = resultFile;
		this.pdbFile = pdbFile;
		parseResult();
		loadStructure();
	}
	
	public IDs2StructureMapper(List<String> IDs1, List<String> IDs2,  String pdbFile){
		this.pdbFile = pdbFile;
		this.IDs1 = IDs1;
		this.IDs2 = IDs2;
		loadStructure();
	}
	
	private void parseResult(){
		this.IDs1 = new ArrayList<String>();
		this.IDs2 = new ArrayList<String>();
		List<String> results = FileIOUtils.createListFromFile(this.resultFile);
		for(int i = 0; i < results.size(); i++){
			String line = results.get(i);
			String[] tokens = line.split("\\s+");
			IDs1.add(tokens[idIndex1]);
			IDs2.add(tokens[idIndex2]);
		}
	}
	
	private void loadStructure(){
		PDBFileReader pdbreader = new PDBFileReader();
        try{
        		//AtomCache cache = new AtomCache();
        		//cache.getFileParsingParams().setAlignSeqRes(true);
        		pdbreader.getFileParsingParameters().setAlignSeqRes(true);	
        		this.structure = pdbreader.getStructure(this.pdbFile);
                System.out.println(structure);
        } catch (IOException e) {
                e.printStackTrace();
        }
        
	}
	
	//map peptide ID to protein structure in files
	private void mapIDs(){
		this.seqs = getChainSeqs();
		this.idMap = new HashMap();
		//mapp the peptide id to pdb chain
		for(int i = 0; i < IDs1.size(); i++){
			String ID1 = this.IDs1.get(i);
			String ID2 = this.IDs2.get(i);
			mapID(ID1);
			mapID(ID2);
		}
		
		for(int i = 0; i < IDs1.size(); i++){
			String ID1 = this.IDs1.get(i);
			String ID2 = this.IDs2.get(i);
			List positions1 = this.idMap.get(ID1);
			List positions2 = this.idMap.get(ID2);
			System.out.print("ID: " + ID1 + "--" + ID2 + " ID1: ");
			for(int j = 0; j < positions1.size()-2; j=j+3){
				System.out.print("("+positions1.get(j) + " :" + positions1.get(j+1) + ", " + positions1.get(j+2) +"); ");
			}
			System.out.println();
			System.out.print("ID: " + ID1 + "--" + ID2 + " ID2: ");
			for(int j = 0; j < positions2.size()-2; j=j+3){
				System.out.print("("+positions2.get(j) + " :" + positions2.get(j+1) + ", " + positions2.get(j+2) +"); ");
			}
			System.out.println();
		}
		
	}
	
	//create external mapping
	private void setIDMap(Map<String, List> idMap){
		this.idMap = idMap;
	}
	
	private void getDistance(){
		//mapp the peptide id to pdb chain
		for(int i = 0; i < IDs1.size(); i++){
			String ID1 = this.IDs1.get(i);
			String ID2 = this.IDs2.get(i);
			if(ID1.compareTo(ID2) > 0){
				String temp = ID1;
				ID1 = ID2;
				ID2 = temp;
			}
			List positions1 = this.idMap.get(ID1);
			List positions2 = this.idMap.get(ID2);
			for(int j = 0; j < positions1.size()-2; j=j+3){
				for(int k = 0; k < positions2.size()-2; k=k+3){
					System.out.print("ID: " + ID1 + "--" + ID2 + "\t: ");
					System.out.print("("+positions1.get(j) + " :" + positions1.get(j+1) + ", " + positions1.get(j+2) +");\tand\t");
					System.out.print("("+positions2.get(k) + " :" + positions2.get(k+1) + ", " + positions2.get(k+2) +"); ");
					System.out.print("\twith distance:"+"\t");
					Chain c1 = this.structure.getChain((Integer)positions1.get(j));
					Chain c2 = this.structure.getChain((Integer)positions2.get(k));
					int index = (Integer)positions1.get(j+2);
					int index2 = (Integer)positions2.get(k+2);
					if(index >= c1.getSeqResLength() || index2 > c2.getSeqResLength()){
						System.out.print("ID: " + ID1 + "--" + ID2 + "\t: ");
						System.out.print("(-1: 0, 0);\tand\t");
						System.out.print("(-1: 0, 0)\t");
						System.out.print("\twith distance:\t1000");
						System.out.print("\n");
						break;
					}
					//System.out.println("seqres: " + c1.getSeqResLength());
					//System.out.println("seqres: " + c2.getSeqResLength());
					Group g1 = c1.getSeqResGroup(index+1);
					Group g2 = c2.getSeqResGroup(index2+1);
					System.out.print("Group1 : " + index + "\t" + g1.getPDBName() + "\t");
					System.out.print("Group2 : " + index2 + "\t" + g2.getPDBName() +"\t");
					try{
						Atom a1 = g1.getAtom("CA");
						Atom a2 = g2.getAtom("CA");
						double[] coor1 = a1.getCoords();
						double[] coor2 = a2.getCoords();
						double dist1 = coor1[0]-coor2[0];
						double dist2 = coor1[1]-coor2[1];
						double dist3 = coor1[2]-coor2[2];
						System.out.print(Math.sqrt(dist1*dist1+dist2*dist2+dist3*dist3));
					}catch(StructureException se){
						
					}
					System.out.println();
				}
			}
			if(positions1.size() == 0 || positions2.size() == 0){
				System.out.print("ID: " + ID1 + "--" + ID2 + "\t: ");
				System.out.print("(-1: 0, 0);\tand\t");
				System.out.print("(-1: 0, 0)\t");
				System.out.print("\twith distance:\t1000");
				System.out.print("\n");
			}
		}
	}
	
	private void mapID(String ID){
		if(ID.startsWith("r")){
			ID = ID.substring(1);
		}
		if(this.idMap.containsKey(ID)){
			return;
		}
		this.idMap.put(ID, new ArrayList<Integer>());
		List<Integer> positions = this.idMap.get(ID);
		for(int j = 0; j < this.seqs.size(); j++){
			String chainSeq = this.seqs.get(j);
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
	
	private List<String> getChainSeqs(){
		List<String> seqs = new ArrayList<String>();
		List<Chain> chainList = this.structure.getChains();
		for(int i = 0; i < chainList.size(); i++){
			Chain chain = chainList.get(i);
			//System.out.println("chain seq in RES: " + chain.getAtomSequence());
			seqs.add(chain.getAtomSequence());
		}
		return seqs;
	}
	
	public List<String> getAllIDs(){
		List<String> IDs = new ArrayList<String>();
		IDs.addAll(this.IDs1);
		IDs.addAll(this.IDs2);
		return IDs;
	}
	
	public static void testID2StrcutMapper(){
		String resultFile = "../mixture_linked/testAnnotation.txt";
		String pdbFile = "../mixture_linked/PDBFiles/1ZAH.pdb";
		IDs2StructureMapper mapper = new IDs2StructureMapper(resultFile, pdbFile);
		mapper.mapIDs();
		mapper.getDistance();
	}
	
	public static void testID2StrcutMapperWithHomolog(){
		String resultFile = "../mixture_linked/testAnnotation.txt";
		String pdbFile = "../mixture_linked/PDBFiles//3UNE.pdb";
		String proteinFasta = "..//mixture_linked//database//Rabbit_uniprot_proteasome.fasta";
		String pdbFasta = "..//mixture_linked//database//3UNE.fasta.txt";
		ProteinIDExtractor protMap = new ProteinIDExtractor(proteinFasta, resultFile);
		HomologMapper hMap = new HomologMapper(proteinFasta, pdbFasta);
		IDs2StructureMapper mapper = new IDs2StructureMapper(resultFile, pdbFile);
		mapper.setIDMap(hMap.createHomologPositionMap(mapper.getAllIDs(), protMap.getPositonMap()));
		mapper.getDistance();
	
	}
	
	public static void main(String[] args){
		//testID2StrcutMapper();
		testID2StrcutMapperWithHomolog();
	}
	
	
}
