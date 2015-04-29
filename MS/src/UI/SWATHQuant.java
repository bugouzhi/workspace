package UI;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.Spectrums.AnnotatedSpectrum;
import dia.SWATHUtils;
import org.Spectrums.TDAStat;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
/**
 * Command line entry points to generate quantitative information from swath data
 * @author Jian
 *
 */
public class SWATHQuant {
	public static void testGenerateIonLibrary(String specLibFile, String resultFile, String fastaFile, String outFile, String outFormat){
//		specLibFile = "../mixture_linked/APSWATH_PPP2RA_MSGFDB_IDs_1pepFDR_plusDecoy2.mgf";
//		resultFile = "../mixture_linked/t0";
//		fastaFile = "../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta";
//		outFile = "..//mixture_linked//testPeakViewLib_out.txt";
		//String ssmResultFile = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPSEcoli_libREp234_msplit_7_pep.txt";
		//getIonLibrary(specLibFile);
		SWATHUtils.getIonLibrary(specLibFile, "", resultFile, fastaFile, outFile, outFormat);
		//getIonLibraryFromSWATH(specLibFile, spectrumFile, resultFile);
	}
	
	public static void testGenerateSpectCount(String resultFile, String fastaFile, String ExpName, String BaitName, String outFile){
//		resultFile = "../mixture_linked/t0";
//		fastaFile = "../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta";
//		BaitName = "PPP2RABait1";
//		ExpName = "PPP2RAExp1";
//		outFile = "../mixture_linked/testSAINTIn.txt";
		TDAStat results = new TDAStat(resultFile, 4, 8, 32, -1);
		results.mapProtein(fastaFile);
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		for(Iterator<String> it = results.proteinMap.keySet().iterator(); it.hasNext();){
			AnnotatedSpectrum protS = results.proteinMap.get(it.next());
			if(((Integer)protS.getAnnotation().get("uniquePepCount")) > 0 ){
				//System.out.println("uniqueCount " + protS.getAnnotation().get("uniquePepCount"));
				out.println(BaitName + "\t" + ExpName + "\t" + protS.protein + "\t" + protS.getAnnotation().get("specCount"));
			}
		}
		out.flush();
		out.close();
	}
	
	public static void generateProhitsOut(String resultFile, String fastaFile, String ExpName, String BaitName, String outFile){
//		resultFile = "../mixture_linked/t0";
//		fastaFile = "../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta";
//		BaitName = "PPP2RABait1";
//		ExpName = "PPP2RAExp1";
//		outFile = "../mixture_linked/testSAINTIn.txt";
		TDAStat results = new TDAStat(resultFile, 4, 8, 7, -1);
		results.mapProtein(fastaFile);
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		List<String> noProteins = new ArrayList<String>();
		noProteins.add("PROTEIN");
		out.println("#Bait\tExperiment\tPeptide\tProtein\tSpecCount\tSSM-Score\n");
		for(Iterator<String> it = results.peptideMap.keySet().iterator(); it.hasNext();){
			AnnotatedSpectrum pepS = results.peptideMap.get(it.next());
			List<String> proteins = (List<String>)pepS.getAnnotation().get("proteins");
			if(proteins == null){
				proteins = noProteins;
			}
			for(int i = 0; i < proteins.size(); i++){
				pepS.score = pepS.score*-1;
				if(pepS.score > 1){
					pepS.score = 0.991;
				}
				out.println(BaitName + "\t" + ExpName + "\t" + pepS.peptide +"\t" + proteins.get(i) +"\t"+pepS.getAnnotation().get("specCount")+"\t"+pepS.score);
			}
		}
		out.flush();
		out.close();
	}
	
	public static Options createOptions(){
		Option help = new Option("help", "help");
		Option resultFile = OptionBuilder.withArgName( " resultFile " )
                .hasArg()
                .withDescription(  "search result file" )
                .isRequired()
                .create( "r" );
		
		Option libFile = OptionBuilder.withArgName( " libFile " )
                .hasArg()
                .withDescription(  "spectral library file" )
                .isRequired()
                .create( "l" );
		
		Option fastaFile = OptionBuilder.withArgName( " fastaFile" )
                .hasArg()
                .withDescription(  "fasta sequence file" )
                .isRequired()
                .create( "d" );
		
		Option outFile = OptionBuilder.withArgName( " output name " )
                .hasArg()
                .withDescription(  "output name/path" )
                .isRequired()
                .create( "o" );
		
		
		Option saintExport = OptionBuilder.withArgName( " 1 | 0 " )
                .hasArg()
                .withDescription(  "protein spectral count with SAINT input format" )
                .create( "SAINTInput" );
		
		Option prohitsExport = OptionBuilder.withArgName( " 1 | 0 " )
                .hasArg()
                .withDescription(  "peptide spectral count with Prohits input format" )
                .create( "ProhitsInput" );
		
		Option expname = OptionBuilder.withArgName( " Exp name " )
                .hasArg()
                .withDescription(  "A name for the experiment (e.g EIF4A2_REP1 or CONTROL_REP2)" )
                .create( "Expname" );
		
		Option baitname = OptionBuilder.withArgName( " bait name " )
                .hasArg()
                .withDescription(  "Name for the bait protein, for user-defined control use CONTROL" )
                .create( "Baitname" );
		
		Option peakviewExport = OptionBuilder.withArgName( " 1 | 0 " )
                .hasArg()
                .withDescription(  "export library in PeakView format" )
                .create( "PeakviewInput" );
		
		Option skylineExport = OptionBuilder.withArgName( " 1 | 0 " )
                .hasArg()
                .withDescription(  "export library in Skyline format" )
                .create( "SkylineInput" );
		
		Option openswathExport = OptionBuilder.withArgName( " 1 | 0 " )
                .hasArg()
                .withDescription(  "export library in Openswath format" )
                .create( "OpenswathInput" );
		
		Options options = new Options();
		options.addOption(help);
		options.addOption(resultFile);
		options.addOption(outFile);
		options.addOption(fastaFile);
		options.addOption(libFile);
		options.addOption(saintExport);
		options.addOption(expname);
		options.addOption(baitname);
		options.addOption(peakviewExport);
		options.addOption(skylineExport);
		options.addOption(openswathExport);
		options.addOption(prohitsExport);
		return options;
	}
	
	
	public static void main(String[] args){
		//legacy format
		if(args[0].equals("PeakViewInput")){
			testGenerateIonLibrary(args[1], args[2], args[3], args[4], "PeakView");
		}else if(args[0].equals("SAINTInput")){
			testGenerateSpectCount(args[1], args[2], args[3], args[4], args[5]);
		}else if(args[0].equals("OpenswathInput")){
			testGenerateIonLibrary(args[1], args[2], args[3], args[4], "OpenSwath");
		}else if(args[0].equals("SkylineInput")){
			testGenerateIonLibrary(args[1], args[2], args[3], args[4], "Skyline");
		}else{
		//new posix compliant commnd line options
		Parser parser = new GnuParser();
	    Options options = createOptions();
	    try {
	    	// parse the command line arguments
	    	CommandLine line = parser.parse( options, args );
	    	if(line.hasOption("help")){
	    			HelpFormatter formatter = new HelpFormatter();
		        	formatter.printHelp( "SWATHQuant", options );
		        	return;
		    }
	    	String resultFile = line.getOptionValue("r");
	    	String libFile = line.getOptionValue("l");
	    	String outFile = line.getOptionValue("o");
	    	String fastaFile = line.getOptionValue("d");
	    	
	    	outFile = Utils.FileIOUtils.generateOutFile(outFile, resultFile, "");
	    	String sep ="_";
	    	
	    	String saintOut = outFile + sep + "SAINT_specCount_out.txt";
	    	//System.out.println("fasta " + fastaFile);
	    	if(line.hasOption("SAINTInput") && line.getOptionValue("SAINTInput").equals("1")){
	    		String expname = "EXP";
	    		String baitname = "BAIT";
	    		if(line.hasOption("Baitname")){
	    			baitname = line.getOptionValue("Baitname");
	    		}
	    		if(line.hasOption("Expname")){
	    			expname = line.getOptionValue("Expname");
	    		}
	    		testGenerateSpectCount(resultFile, fastaFile, expname, baitname, saintOut);
	    	}
	    	
	    	String prohitsOut = outFile + sep + "prohits_specCount_out.txt";
	    	if(line.hasOption("ProhitsInput") && line.getOptionValue("ProhitsInput").equals("1")){
	    		String expname = "EXP";
	    		String baitname = "BAIT";
	    		if(line.hasOption("Baitname")){
	    			baitname = line.getOptionValue("Baitname");
	    		}
	    		if(line.hasOption("Expname")){
	    			expname = line.getOptionValue("Expname");
	    		}
	    		generateProhitsOut(resultFile, fastaFile, expname, baitname, prohitsOut);
	    	}
	    	
	    	String peakviewOut = outFile + sep + "Peakview_assaylib.txt";
	    	String skylineOut = outFile + sep + "Skyline_assaylib.txt";
	    	String openswathOut = outFile + sep + "Openswath_assaylib.txt";    	
	    	if(line.hasOption("PeakviewInput") && line.getOptionValue("PeakviewInput").equals("1")){
	    		testGenerateIonLibrary(libFile, resultFile, fastaFile, peakviewOut, "PeakView");
	    	}
	    	if(line.hasOption("SkylineInput") && line.getOptionValue("SkylineInput").equals("1")){
	    		testGenerateIonLibrary(libFile, resultFile, fastaFile, skylineOut, "Skyline");
	    	}
	    	if(line.hasOption("OpenswathInput") && line.getOptionValue("OpenswathInput").equals("1")){
	    		testGenerateIonLibrary(libFile, resultFile, fastaFile, openswathOut, "OpenSwath");
	    	}
		 }
		   catch( ParseException exp ) {
		        // oops, something went wrong
		        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
		        HelpFormatter formatter = new HelpFormatter();
	        	formatter.printHelp( "SWATHQuant", options);
		   }
		
		}
	}
}
