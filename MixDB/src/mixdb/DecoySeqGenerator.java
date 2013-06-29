package mixdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class DecoySeqGenerator {
	private String fastaFile;
	private String outFile;
	public DecoySeqGenerator(String fastaFile){
		this.fastaFile = fastaFile;
		this.outFile = this.fastaFile.substring(0, this.fastaFile.lastIndexOf('.')) 
				+ "_plusDecoy.fasta";
		System.out.println("out: " + outFile);
	}
	
	public DecoySeqGenerator(String fastaFile, String outFile){
		this.fastaFile = fastaFile;
		this.outFile = outFile;
		System.out.println("out: " + outFile);
	}
	
	public void generateDecoy(){
		try{
			BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
			String line = reader.readLine();
			String header = "";
			String seq ="";
			while(line != null){
				//System.out.println("line is " +  line);
				if(line.startsWith(">")){
					if(seq.length() > 0){
						//System.out.println("writing" + header);
						writer.write(header+"\n");
						writer.write(seq + "\n");
						writer.write(">X_" + header.substring(1) + "\n");
						writer.write(new StringBuffer(seq).reverse().toString() + "\n");
					}
					header = line;
					seq = "";
				}else{
					if(line.length() > 0){
						seq = seq + line;
					}
				}
				line = reader.readLine();
			}
			if(seq.length() > 0){
				//System.out.println("writing" + header);
				writer.write(header+"\n");
				writer.write(seq + "\n");
				writer.write(">X_" + header.substring(1) + "\n");
				writer.write(new StringBuffer(seq).reverse().toString() + "\n");
			}
			reader.close();
			writer.flush();
			writer.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public static void main(String[] args){
		if(args.length != 2){
			System.out.println("usage: java -Xmx1000M -jar DecoySeqGenerate.jar <fasta File> <output file>");
			return;
		}
		try{
			String fasta = args[0];
			String out = args[1];
			//String fasta = "../mixture_linked/Ecoli_genome.fasta";
			DecoySeqGenerator decoy = new DecoySeqGenerator(fasta, out);
			decoy.generateDecoy();
		}catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
