package org.Spectrums;
/**
 * 
 * Given two set of sequence we try to map sequence in one
 * set to sequence in another set. Mapping is very simple,
 * it is the highest-scoring matches between the two sets.
 * Score compute by alignment score.
 * @author Jian Wang
 *
 */
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SimpleSubstitutionMatrix;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;

public class HomologMapper {
	private String seqFile1;
	private String seqFile2;
	private LinkedHashMap<String, ProteinSequence> seq1;
	private LinkedHashMap<String, ProteinSequence> seq2;
	private LinkedHashMap<String, List<AlignedSequence>> homologMap;
	private Map<String, Integer> chainMap;
	
	public HomologMapper(String seqFile1, String seqFile2){
		this.seqFile1 = seqFile1;
		this.seqFile2 = seqFile2;
		try{	
			this.seq1 = FastaReaderHelper.readFastaProteinSequence(new File(this.seqFile1));
			this.seq2 = FastaReaderHelper.readFastaProteinSequence(new File(this.seqFile2));
		}catch(Exception ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		getChainIndex();
		computeAlignment();
	}
	
	private void getChainIndex(){
		this.chainMap = new HashMap<String, Integer>();
		int chainIndex = 0;
		for(Iterator<ProteinSequence> it = seq2.values().iterator(); it.hasNext();){
			ProteinSequence prot = it.next();
			System.out.println("homolog name: " + prot.getOriginalHeader());
			chainMap.put(prot.getOriginalHeader(), chainIndex);
			chainIndex++;
		}
	}
	
	private void computeAlignment(){
		this.homologMap = new LinkedHashMap<String, List<AlignedSequence>>();
        SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
        GapPenalty gapScores = new SimpleGapPenalty();
        List<ProteinSequence> seqList = new ArrayList();
        List<List<AlignedSequence>> bestPair = new ArrayList();
        for(Iterator<ProteinSequence> it = seq1.values().iterator();  it.hasNext();){
        	ProteinSequence s1 = it.next();
        	ProteinSequence bestMatch = null;
        	Profile bestAlign = null;
      		int bestScore = -1000000;
       		for(Iterator<ProteinSequence> it2 = seq2.values().iterator(); it2.hasNext();){
       			ProteinSequence s2 = it2.next();       			
       			PairwiseSequenceAligner align = Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.GLOBAL, gapScores, matrix);
       			int currentScore = align.getScore();
       			bestMatch = currentScore > bestScore ? s2 : bestMatch;
       			bestAlign = currentScore > bestScore ? align.getProfile() : bestAlign;
       			bestScore = currentScore > bestScore ? currentScore : bestScore;
       			//break;
       		}
       		System.out.println(s1.getOriginalHeader() + "\tbest score:\t" + bestScore +"\t@seq: " + bestMatch.getOriginalHeader());
       		System.out.println("Aligned sequence: " + bestAlign.getAlignedSequences().size());
       		for(int i = 0; i < bestAlign.getAlignedSequences().size(); i++){
       			System.out.println(bestAlign.getAlignedSequences().get(i));
       		}
       		//AlignedSequence orig = (AlignedSequence)bestAlign.getAlignedSequences().get(0);
       		//System.out.println("original: " + ((ProteinSequence)orig.getOriginalSequence()).getOriginalHeader());
       		System.out.println();
       		System.out.println("putting key: " + s1.getOriginalHeader() + "\t" + s1.getOriginalHeader().length());
       		if(bestScore > 20){
       			this.homologMap.put(s1.getOriginalHeader(), bestAlign.getAlignedSequences());
       		}
        }
        
	}
	
	private String getChainID(String name){
		return name.split("|")[0].split(":")[1];
	}
	
	public String getEquivalentSeq(String proteinName, String peptide){
		ProteinSequence prot = null;
		if(!this.seq1.containsKey(proteinName)){
			System.err.println("warning: cannot find protein sequences");
		}else{
			prot = this.seq1.get(proteinName);
		}
		String protStr = prot.getSequenceAsString();
		int startIndex = protStr.indexOf(peptide);
		if(startIndex < 0){
			System.err.println("warning: peptide cannot be mapped to sequence");
			return "";
		}
		AlignedSequence align = this.homologMap.get(proteinName).get(0);
		AlignedSequence align2 = this.homologMap.get(proteinName).get(1);
		return align2.getSequenceAsString().substring(gappedIndex(align, startIndex), 
					gappedIndex(align, startIndex+peptide.length()));
		
	}
	
	
	public Object[] getEquivalentPosition(String proteinName, int index){
		Object[] position = new Object[2];
		ProteinSequence prot = null;
		if(!this.homologMap.containsKey(proteinName)){
			System.err.println("warning: cannot find protein sequences");
			return position;
		}
		AlignedSequence align = this.homologMap.get(proteinName).get(0);
		AlignedSequence align2 = this.homologMap.get(proteinName).get(1);
		String homologName = ((ProteinSequence)align2.getOriginalSequence()).getOriginalHeader();
		int gapInd = gappedIndex(align, index);
		int nonGap = nonGappedIndex(align2, gapInd);
		System.out.println("homoglog name: " + homologName);
		System.out.println("index: " + index + "\t" + gapInd + "\t" + nonGap);
		position[0] = this.chainMap.get(homologName);
		position[1] = nonGap;
		return position;
	}
	
	//given non-gapped index, return gapped index in alignment
	private int gappedIndex(AlignedSequence seq, int index){
		int count = 0;
		int i = 0;
		for(; i < seq.getLength(); i++){
			if(!seq.isGap(i)){
				count++;
			}
			if(count >= index){
				return i;
			}
		}
		return -1;
	}
	
	//given gapped index, return first non-gapped index in alignment    //should we consider first gapped or cloest non-gapped?
	private int nonGappedIndex(AlignedSequence seq, int index){
		int count = 0;
		int i = 0;
		for(; i < index; i++){
			if(!seq.isGap(i)){
				count++;
			}
		}
		return count;
	}

	
	public Map<String, List> createHomologPositionMap(List<String> ids, Map<String, List<Object>> positionMap){
		Map<String, List> hIdMap = new HashMap<String, List>();
		System.out.println("position map size: " + positionMap.keySet().size());
		for(Iterator<String> idIter = positionMap.keySet().iterator(); idIter.hasNext();){
			String id = idIter.next();
			System.out.println("ID is: " + id);
			List<Object> position = positionMap.get(id);
			List<Object> homologPos = new ArrayList<Object>();
			hIdMap.put(id, homologPos);
			//System.out.println("position List size: " + position.size());
			for(int i = 0; i < position.size()-1; i+=2){
				System.out.println(position.get(i));
				System.out.println(position.get(i+1));
				Object[] hPosition = this.getEquivalentPosition((String)position.get(i), ((Integer)position.get(i+1)).intValue());
				if(hPosition[0] != null){
					homologPos.add(hPosition[0]);
					homologPos.add(hPosition[1]);
					homologPos.add((Integer)hPosition[1]+id.indexOf('K'));
					System.out.println("mapping: " + id + " native: " + position.get(i) + ", " + position.get(i+1) + "\tmapped to\t:" + hPosition[0] + ", " + hPosition[1]);
					System.out.println("native protein length: " + ((String)position.get(i)).length());
				}
			}
		}
		return hIdMap;
	}
	
	public static void testHomologMapper(){
		String file1 = "..//mixture_linked//database//Rabbit_uniprot_proteasome.fasta";
		String file2 = "..//mixture_linked//database//3UNE.fasta.txt";
		HomologMapper map = new HomologMapper(file1, file2);
	}
	
	public static void main(String[] args){
		testHomologMapper();
	}
	
}
