package IO;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;

import org.Spectrums.AnnotatedSpectrum;
import org.Spectrums.TDAStat;

/**
 * Generate peptide prophet output so can be used to run protein prophet
 * 
 * @author Jian
 *
 */


public class PeptideProphetGenerator {
	String headerFile="../mixture_linked/PeptideProphetHeader.xml";
	String resultFile;
	String outFile;
	String fastaFile;
	PrintStream out = Utils.FileIOUtils.getOutStream("");
	
	public PeptideProphetGenerator(String resultFile, String fastaFile){
		this.resultFile = resultFile;
		File f = new File(fastaFile);
		if(f.exists()){
			this.fastaFile = f.getAbsolutePath();
		}else{
			throw new IllegalArgumentException("Fasta file does ot exist");
		}
	}
	
	public void generateProphetFile(String outFile){
		this.outFile = outFile;
		this.out = Utils.FileIOUtils.getOutStream(outFile);
		this.printHeader();
		this.printPeptideResult();
		this.printCloseStatement();
		out.flush();
		out.close();
	}
	
	public void printHeader(){
		List<String> headers = Utils.FileIOUtils.createListFromFile(this.headerFile);
		for(int i = 0; i < headers.size(); i++){
			out.println(headers.get(i));
		}
		out.println("<analysis_timestamp analysis=\"peptideprophet\" time=\"2014-08-25T16:14:09\" id=\"1\"/>");
		out.println("<analysis_timestamp analysis=\"database_refresh\" time=\"2014-08-25T16:14:09\" id=\"1\">");
		out.println("<database_refresh_timestamp database=\"" + this.fastaFile +  "\" min_num_enz_term=\"1\"/>");
		out.println("</analysis_timestamp>");
	}
	
	public void printPeptideResult(){
		TDAStat stat = new TDAStat(resultFile, 1, -1, 7, 6, 8, 11, 1, false, 0.01, this.fastaFile);
		int count = 0;
		for(Iterator<AnnotatedSpectrum> it = stat.peptideMap.values().iterator(); it.hasNext();){
			 AnnotatedSpectrum s = it.next();
			 if((Double)s.getAnnotation().get("pepfdr") < 0.05){
				 printPeptideMatch(s, count);
			 	count++;
			 	if(count > 1000000){
			 		break;
			 	}
			 }
		}
	}
	
	public void printPeptideMatch(AnnotatedSpectrum s, int counter){
		List<String> proteins = (List<String>)s.getAnnotation().get("proteins");
		double score = s.score*100000000;//(Double)s.getAnnotation().get("pepfdr") + Math.random();
		System.out.println("score is: " + s.score + "\t" + score);
		score = score > 1 ? 1.0 : score;
		if(proteins == null || proteins.size() == 0){
			return;
		}
		String protein;
		if(counter > 5 && counter < 10) 
			protein = "gi|90111377|ref|NP_416564.4|";
		else 
			protein = proteins.get(0).split("\\s+")[0];
		out.print("<spectrum_query spectrum=\"SpectrumFile.00000.00000." + s.scanNumber +  "\"" 
		        +" start_scan=\"" + s.scanNumber  + "\"" 
				+ " end_scan=\"" + s.scanNumber + "\"" 
				+ " precursor_neutral_mass=\"" + s.parentMass + "\"" 
				+ " assumed_charge=\"" + s.charge + "\"" 
				+ " index=\"" + counter
				+ "\" retention_time_sec=\"1000.0\">" + "\n" 
		+ "<search_result>\n" 

		+ "<search_hit hit_rank=\"1\" peptide=\"" + s.peptide + "\""  
		+ " peptide_prev_aa=\"R\" peptide_next_aa=\"G\" protein=\"" + protein/*proteins.get(0).split("\\s+")[0]*/ + "\"" 
		+ " num_tot_proteins=\"" + proteins.size() + "\" num_matched_ions=\"20\" tot_num_ions=\"20\" calc_neutral_pep_mass=\"" + s.parentMass + "\""
		+ " massdiff=\"" + (s.parentMass - s.parentMass) +"\""
		+ " num_tol_term=\"2\" num_missed_cleavages=\"0\" is_rejected=\"0\">\n");
		
		int size = proteins.size();
		for(int i = 1; i < size; i++){
			String alternateprot = proteins.get(i).split("\\s+")[0];
			out.println("<alternative_protein protein=\"" + alternateprot + "\"");
		}
		double prob1 = 0.0;
		double prob2 = Math.random();
		double prob3 = 1-score;
		out.println("<analysis_result analysis=\"peptideprophet\">");
		out.println("<peptideprophet_result probability=\"" + (1-score) + "\"" 
				     + " all_ntt_prob=\"(" + prob1 + ", " + prob2 + ", " + prob3 +  ")\">");
		out.print("</peptideprophet_result>\n" 
		           + "</analysis_result>\n"
		           + "</search_hit>\n"
		           + "</search_result>\n"
		           + "</spectrum_query>\n\n\n");
	}
	
	public void printCloseStatement(){
		out.println("</msms_run_summary>");
		out.println("</msms_pipeline_analysis>");
	}
	
	public static void testGenerateProphetOutFile(){
		String resultFile = "../mixture_linked/SWATH/Swath_Human_searches/IDAMsgfdb/UPS_Human_IDAREP3_msgfdb_1.txt";
		String fastaFile = "../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta";
		String outFile = "../mixture_linked/testProphetout.xml";
		PeptideProphetGenerator gen = new PeptideProphetGenerator(resultFile, fastaFile);
		gen.generateProphetFile(outFile);
	}
	
	public static void main(String[] args){
		testGenerateProphetOutFile();
	}
}	



