package org.Spectrums;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import org.biojava.bio.proteomics.MassCalc;
import org.biojava.bio.proteomics.Digest;
import org.biojava.bio.proteomics.Protease;
import org.biojava.bio.proteomics.ProteaseManager;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.BioException;
import org.biojava.utils.ChangeVetoException;

/**
 * Digest a proteins into peptides
 * @author jian wang
 *
 */
public class ProteinDigester {
	public static Digest digester = new Digest();
	public static List<String> digestProteinByTrypsin(String protein, String name){
		return digestProteinByTrypsin(protein, name, 6, 25);
	}
	
	public static List<String> digestProteinByTrypsin(String protein, String name, int minLength, int maxLength){
		List<String> peptides = new ArrayList();
		try{
			digester.setProtease(ProteaseManager.getProteaseByName(Protease.CHYMOTRYP));
			Sequence p = ProteinTools.createProteinSequence(protein, name);
			digester.setSequence(p);
			digester.setMaxMissedCleavages(500);
			digester.addDigestFeatures();
			Iterator it = p.features();
			Feature f;
			String pep;
			while(it.hasNext()){
				f = (Feature)it.next();
				//System.out.println(f.getSymbols().seqString());
				if(f.getSymbols().length() >= minLength && f.getSymbols().length() <= maxLength){
					System.out.println(f.getSymbols().seqString());
					//peptides.add(f.getSymbols().seqString());
				}
			}
			//mass = mcalc.getMolecularWeight(ProteinTools.createProtein(pep));
		}catch(IllegalSymbolException e){
			System.out.println("Peptide not valid: " + protein);
			System.out.println(e.getMessage());
		}catch(BioException e){
			System.out.println(e.getMessage());
		}catch(ChangeVetoException e){
			System.out.println(e.getMessage());
		}
		return peptides;
	}
	
	public static void digestProteins(String fasta){
		try{
			BufferedReader bf = new BufferedReader(new FileReader(fasta));
			String currentLine = bf.readLine();
			StringBuffer protein = null;
			String header = null;
			while(currentLine != null){
				//System.out.println("line is: " + currentLine);
				if(currentLine.startsWith(">")){
					if(protein != null){
						ProteinDigester.digestProteinByTrypsin(protein.toString(), header);
						//System.out.println(header);
						//System.out.println(protein.toString());
						
						protein = new StringBuffer();
					}
					header = currentLine;	
					protein = new StringBuffer();
				}else{
					protein.append(currentLine);
				}
				currentLine = bf.readLine();
			}
			if(protein != null){
				ProteinDigester.digestProteinByTrypsin(protein.toString(), header);
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void testDigester(String protein){
		//digestProteinByTrypsin(protein, "name");
		//digestProteins("..\\mixture_linked\\tryp_unlinked_proteins.fasta");
		//digestProteins("..\\mixture_linked\\Ecoli_genome.fasta");
		//digestProteins("..\\mixture_linked\\Ecoli_genomeReverse.fasta");
		digestProteins("..\\mixture_linked\\database\\Human_allproteins.fasta");
	}
	
	
	public static void main(String[] args){
		String protein = "MQVSVETTQGLGRRVTITIAADSIETAVKSELVNVAKKVRIDGFRKGKVPMNIVAQRYGASVRQDVLGDL" 
			+ "MSRNFIDAIIKEKINPAGAPTYVPGEYKLGEDFTYSVEFEVYPEVELQGLEAIEVEKPIVEVTDADVDGM" 
			+ "LDTLRKQQATWKEKDGAVEAEDRVTIDFTGSVDGEEFEGGKASDFVLAMGQGRMIPGFEDGIKGHKAGEE"
			+ "FTIDVTFPEEYHAENLKGKAAKFAINLKKVEERELPELTAEFIKRFGVEDGSVEGLRAEVRKNMERELKS"
			+ "AIRNRVKSQAIEGLVKANDIDVPAALIDSEIDVLRRQAAQRFGGNEKQALELPRELFEEQAKRRVVVGLL"
			+ "LGEVIRTNELKADEERVKGLIEEMASAYEDPKEVIEFYSKNKELMDNMRNVALEEQAVEAVLAKAKVTEK"
			+ "ETTFNELMNQQA";
		testDigester(protein);
	}
	
}
