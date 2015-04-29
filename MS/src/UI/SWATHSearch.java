package UI;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

import dia.SWATHMSPLITSearch;

public class SWATHSearch {
	public static void SWATHSearch(){
		
	}
	
	public static Options createOptions(){
		Option help = new Option("help", "help");
		Option resultFile = OptionBuilder.withArgName( "spectrum file" )
                .hasArg()
                .withDescription(  "search result file" )
                .isRequired()
                .create( "s" );
		
		Option libFile = OptionBuilder.withArgName( "library file" )
                .hasArg()
                .withDescription(  "spectral library file" )
                .isRequired()
                .create( "l" );
		
		Option outFile = OptionBuilder.withArgName( "output file" )
                .hasArg()
                .withDescription(  "output file" )
                .isRequired()
                .create( "o" );
		
		Option ptmFile = OptionBuilder.withArgName( "ptm file" )
                .hasArg()
                .withDescription(  "ptm file" )
                .isRequired()
                .create( "ptm" );
		
		
		Option cycle = OptionBuilder.withArgName( "SWATH cycle" )
                .hasArg()
                .withDescription(  "Number of scans to cover the whole m/z range in the experiment" )
                .create( "Cycle" );
		
		Option isolation = OptionBuilder.withArgName("Isolation Window")
                .withDescription(  "Precusor isolation window in Da" )
                .create( "Cycle" );
		
		Option fragment_tol = OptionBuilder.withArgName("mass tolerance")
                .withDescription(  "Fragment mass tolerance in ppm" )
                .create( "Cycle" );
		
		Options options = new Options();
		options.addOption(help);
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
	        	formatter.printHelp( "SWATHSearch", options );
	        	return;
	        } 
	        SWATHMSPLITSearch search = new SWATHMSPLITSearch();
			if(line.hasOption("precursor"))
				search.parent = Double.parseDouble(line.getOptionValue("precursor"));
			if(line.hasOption("fragment"))
				search.fragment = Double.parseDouble(line.getOptionValue("fragment"));
			search.queryFile = line.getOptionValue("s");
			if(line.hasOption("cycle"))
					search.SWATHCycle = Integer.parseInt(line.getOptionValue("cycle"));
			search.libraryFile = line.getOptionValue("l");
			search.outFile = line.getOptionValue("o");

	    }
	    catch( ParseException exp ) {
	        // oops, something went wrong
	        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
	        HelpFormatter formatter = new HelpFormatter();
        	formatter.printHelp( "SWATHFilter", options);
	    }
	}
}
