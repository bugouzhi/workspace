package UI;
/**
 * Map peptide to a given proteins
 * @author Jian
 *
 */

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.Spectrums.ProteinIDExtractor;

public class ProteinMapper {
	int pepInd;
	int protInd;
	String resultFile;
	ProteinIDExtractor proteinIds;
	public ProteinMapper(String resultFile, String fastaFile, int pepInd, int protInd){
		this.pepInd = pepInd;
		this.protInd = protInd;
		this.resultFile = resultFile;
		proteinIds = new ProteinIDExtractor(fastaFile, resultFile, pepInd);
		
	}
	
	public void reMapProteins(String outFile){
		try{
			Map<String, List<String>> peptideMap = this.proteinIds.getPeptideMap();
			BufferedReader reader = Utils.FileIOUtils.createReaderFromFile(resultFile);
			PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
			String line = reader.readLine();
			while(line != null){
				String[] tokens = line.split("\\t");
				if(line.startsWith("#") || line.startsWith("SpectrumFile")){
					out.println(line + "\tIsUnique\tProteins");
					line = reader.readLine();
					continue;
				}
				String[] peptides = tokens[pepInd].split("!");
				String[] NISTproteins = tokens[protInd].split("!");
				String prot = "";
				String isUnique =  "";
				String prots = "";
				for(int i = 0; i < peptides.length; i++){
					String pep = Utils.StringUtils.getStrippedSeq(peptides[i]);
					List<String> proteins = peptideMap.get(pep);					
					if(peptides.length > 1){
						//if(proteins != null)
							//System.out.println(pep + "\t" + proteins.size());
					}
					if(i > 0){
						//System.out.println("Peptide is: " + pep);
						prot = prot+"!";
						prots = prots + "!";
						isUnique = isUnique +"!";
					}
					if(proteins == null || proteins.size() == 0){
						prot = prot + NISTproteins[i] + "(NIST)";//prot.replaceAll("!", replacement)//prot + tokens[protInd]+"(NIST)";//"NOT-PRESENT";
						prots = prots + NISTproteins[i] + "(NIST)";//prots+tokens[protInd]+"(NIST);";
					}else{
						prot = prot+proteins.get(0).split("\\s+")[0];
						for(int k = 0; k < proteins.size(); k++){
							if(k > 0){
								prots = prots + ";";
							}
							prots=prots+proteins.get(k).split("\\s+")[0];	
						}
					}
					if(proteins == null || proteins.size() > 1){
						isUnique = isUnique + "0";
					}else{
						isUnique = isUnique + "1";
					}
				}
				for(int i = 0; i < tokens.length; i++){
					if(i == protInd){
						out.print(prot + "\t");
					}else{
						out.print(tokens[i] + "\t");
					}
				}
				out.print(isUnique+"\t");
				out.println(prots);
				line = reader.readLine();
			}
			
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	
	private String insertNISTProteins(String protein){
		String[] proteins = protein.split("!");
		String NISTProtein = "";
		for(int i = 0; i < proteins.length; i++){
			proteins[i] = proteins[i] + "NIST";
		}
		
		for(int i = 0; i < proteins.length; i++){
			if(i > 0){
				NISTProtein = NISTProtein + "!";
			}
			NISTProtein = NISTProtein + proteins[i];
		}
		return NISTProtein;
		
	}
	public static void main(String[] args){
		CommandLineParser parser = new CommandLineParser(args);
		String resultFile = parser.getString(0);
		String fastaFile = parser.getString(1);
		int pepInd = parser.getInteger(2);
		int protInd = parser.getInteger(3);
		String outFile = parser.getString(4);
		//resultFile = "..//mixture_linked/Genentech_speclib_search/temp";
		//fastaFile = "..//mixture_linked/database/Human_allproteins.fasta";
		//pepInd = 0;
		//protInd = 1;
		//outFile = "..//mixture_linked//prot_remapp.txt";
		ProteinMapper mapper = new ProteinMapper(resultFile, fastaFile, pepInd, protInd);
		mapper.reMapProteins(outFile);
	}
}
