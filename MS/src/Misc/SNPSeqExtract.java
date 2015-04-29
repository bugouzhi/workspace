package Misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.Spectrums.ProteinIDExtractor;


/**
 * This script extract the set of corresponding SNP seqs from a of identified proteins
 * from a set of PSMs
 * @author Jian
 *
 */
public class SNPSeqExtract {

	/**
	 * This method extract SNP sequences based on identified nonSNP sequence
	 * @param resultFile
	 * @param fastaFile - the database use to identify proteins
	 * @param SNPFile - the SNP database file, the naming convention must be the original protein name in fastaFile plus dbSNP tag
	 * @param protInd
	 * @param out
	 */
	public static void getSNPSeqs(String resultFile, String fastaFile, String SNPFile, int pepInd, String outFile){
		try{
			Set<String> peps = new HashSet<String>();
			BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(resultFile);
			String line = reader.readLine();
			PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
			boolean startPrint = false;
			while(line != null){
				String[] tokens = line.split("\\t");
				peps.add(Utils.StringUtils.stripNeighborRes(tokens[pepInd]));
				line = reader.readLine();
			}
			System.out.println("Identified peptides " + peps.size());
			ProteinIDExtractor protIDs = new ProteinIDExtractor(peps, fastaFile);
			Set<String> prots = protIDs.getProteinMap().keySet();
			System.out.println("Identified proteins: " + prots.size());
//			for(Iterator<String> it=prots.iterator(); it.hasNext();){
//				System.out.println(it.next());
//			}
			Map<String, List<String>> pepMap = protIDs.getPeptideMap();
			int i = 0;
			for(Iterator<String> it=pepMap.keySet().iterator(); it.hasNext();){
				i++;
				String pep = it.next();
				List<String> proteins = pepMap.get(pep);
				if(proteins==null || proteins.size()==0){
					System.out.println("Unmapped peptides: " + pep);
				}
			}
			BufferedReader reader2 = Utils.FileIOUtils.createReaderFromFile(SNPFile);
			line = reader2.readLine();
			while(line != null){
				if(line.startsWith(">")){
					String[] tokens = line.split("-dbSNP");
					String prot = tokens[0].substring(1);
					//System.out.println(prot);
					if(prots.contains(prot)){
						out.println(line);
						startPrint = true;
					}else{
						startPrint = false;
					}
				}else{
					if(startPrint){
						out.println(line);
					}
				}
				line=reader2.readLine();
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	/**
	 * Extract SNP sequences from idetnified proteins
	 * @param protIDsFile
	 * @param SNPFile
	 * @param pepInd
	 * @param outFile
	 */
	public static void getSNPSeqs(String protIDsFile, String SNPFile, String outFile){
		try{
			Set<String> peps = new HashSet<String>();
			BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(protIDsFile);
			String line = reader.readLine();
			PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
			Set<String> prots = new HashSet<String>();
			boolean startPrint = false;
			while(line != null){
				if(line.startsWith(">")){
					String prot = line.substring(1);
					if(prot.startsWith("REV_")){
						prot = prot.substring(4);
					}
					prots.add(prot);
				}
				line = reader.readLine();
			}
			//System.out.println("Identified peptides " + peps.size());
			System.out.println("Identified proteins: " + prots.size());
//			for(Iterator<String> it=prots.iterator(); it.hasNext();){
//				System.out.println(it.next());
//			}

			BufferedReader reader2 = Utils.FileIOUtils.createReaderFromFile(SNPFile);
			line = reader2.readLine();
			while(line != null){
				if(line.startsWith(">")){
					String[] tokens = line.split("-dbSNP");
					String prot = tokens[0].substring(1);
					//System.out.println(prot);
					if(prots.contains(prot)){
						out.println(line);
						startPrint = true;
					}else{
						startPrint = false;
					}
				}else{
					if(startPrint){
						out.println(line);
					}
				}
				line=reader2.readLine();
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	
	public static void main(String[] args){
//		args = new String[5];
//		args[0] = "..//mixture_linked//SunJun_Privacy_Data//out2.txt";
//		args[1] = "..//mixture_linked//SunJun_Privacy_Data//out1.fasta";
//		args[2] = "..//mixture_linked//SunJun_Privacy_Data//out2.fasta";
//		args[3] = "7";
//		args[4] = "..//mixture_linked//SunJun_Privacy_Data//Extracted_SNPs.fasta";
//		getSNPSeqs(args[0], args[1], args[2], Integer.parseInt(args[3]), args[4]);
		getSNPSeqs(args[0], args[1], args[2]);
	}
	
	
	
}
