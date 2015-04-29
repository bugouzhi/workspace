package UI;
import org.Spectrums.TDAStat;
import org.apache.commons.cli.*;

import rt.RTAligner;


/**
 * Provide command line access to filtering of SWATH-MSPLIT/MSPLTI-DIA results
 * @author Jian
 *
 */
public class SWATHFilter{
	public static void filterSWATH(String resultFile, String outFile, String fdr){
		//resultFile = "../mixture_linked/SWATH/testRTConstraints/Human_500nglysate_QTOF5600_tppSpstconsensus_lib_swathmsplit_out_rtFiltered.txt";
		//outFile = "../mixture_linked/SWATH/testRTConstraints/filtered.txt";
		//fdr = "0.01";
		TDAStat.runTDA(resultFile, outFile, Double.parseDouble(fdr), 4, 6, 8, 32 , -1);	
	}
	
	public static void filterSWATHWithRT(String resultFile, String outFile, String fdr, double minRT, double maxRT){
		//resultFile = "..//mixture_linked//SWATH/testRTConstraints//APSWATH_Swathmsplit_withHumanlysateLib_msplit_0_pep.txt";
		//outFile = "..//mixture_linked//SWATH/testRTConstraints//APSWATH_Swathmsplit_withHumanlysateLib_msplit_0_out.txt";
		//fdr = "0.01";
		RTAligner.filterMSPLITDIAWithRT(resultFile, outFile, Double.parseDouble(fdr), minRT, maxRT);
	}
	
	public static Options createOptions(){
		Option help = new Option("help", "help");
		Option resultFile = OptionBuilder.withArgName( "resultFile" )
                .hasArg()
                .withDescription(  "search result file" )
                .isRequired()
                .create( "r" );
		
		Option useRT = OptionBuilder.withArgName( " 0 | 1 " )
                .hasArg()
                .withDescription(  "use Reteintion time to filter result" )
                .create( "rt" );
		
		Option minRT = OptionBuilder.withArgName( "minRT" )
                .hasArg()
                .withDescription(  "min RT used to build RT correlation" )
                .create( "minRT" );
		
		Option maxRT = OptionBuilder.withArgName( "maxRT" )
                .hasArg()
                .withDescription(  "max RT used in build RT corrlation" )
                .create( "maxRT" );
		
		Option outFile = OptionBuilder.withArgName( "filtered output file" )
                .hasArg()
                .withDescription(  "output for filtered result" )
                .isRequired()
                .create( "o" );
		
		Option fdr = OptionBuilder.withArgName( "fdr" )
                .hasArg()
                .withDescription(  "fdr for filtering result" )
                .withType(new Double(1.0))
                .isRequired()
                .create( "fdr" );
		
		Options options = new Options();
		options.addOption(help);
		options.addOption(useRT);
		options.addOption(resultFile);
		options.addOption(outFile);
		options.addOption(fdr);
		options.addOption(minRT);
		options.addOption(maxRT);
		return options;
	}
	
	public static void main( String[] args ) {
	    // create the parser
	    Parser parser = new GnuParser();
    	Options options = createOptions();
	    try {
	        // parse the command line arguments
	        CommandLine line = parser.parse( options, args );
	        if(line.hasOption("help")){
	        	HelpFormatter formatter = new HelpFormatter();
	        	formatter.printHelp( "SWATHFilter", options );
	        	return;
	        }
	        String resultFile = line.getOptionValue("r");
	        String outFile = line.getOptionValue("o");
	        String fdr = line.getOptionValue("fdr");
	        outFile = Utils.FileIOUtils.generateOutFile(outFile, resultFile, "_filtered.txt");
	        try{
	        	Double.parseDouble(fdr);
	        }catch(NumberFormatException e){
	        	System.err.println("fdr is not a valid number");
	        	System.err.println(e.getMessage());
	        	e.printStackTrace();
	        	return;
	        }
	        if(line.hasOption("rt") && line.getOptionValue("rt").equals("1")){
	        	double minRT = 0;
	        	double maxRT = Double.MAX_VALUE;
	        	if(line.hasOption("minRT")){
	        		minRT = Double.parseDouble(line.getOptionValue("minRT"));
	        	}
	        	if(line.hasOption("maxRT")){
	        		maxRT = Double.parseDouble(line.getOptionValue("maxRT"));
	        	}
	        	filterSWATHWithRT(resultFile, outFile, fdr, minRT, maxRT);    	
	        }else{
	        	filterSWATH(resultFile, outFile, fdr);
	        }
	    }
	    catch( ParseException exp ) {
	        // oops, something went wrong
	        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
	        HelpFormatter formatter = new HelpFormatter();
        	formatter.printHelp( "SWATHFilter", options);
	    }
	}

}
