package mixdb;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.Spectrums.LazyEvaluatedSpectrum;
import org.Spectrums.MZXMLReader;
import org.Spectrums.Mass;
import org.Spectrums.MixturePeakScoreLearner;
import org.Spectrums.PTM;
import org.Spectrums.Peptide;
import org.Spectrums.PeptideLite;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.SimpleProbabilisticScorer;
import org.Spectrums.SortedMZXMLReader;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumComparator;
import org.Spectrums.SpectrumLibSearcher;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;

import Utils.FileIOUtils;
/**
 * This class is the main entry point into mixdb tool, it setup input parameters from command line
 * perform different type of searches depending on input and generates outputs.  Note MixDBSearch can
 * perform several different kind of searches.  However except for the standard searchDB() methods, the
 * others are consider experimental methods in this stage.  As they have not be tested on many datasets,
 * and thus may need further test/validations for them in future development
 * @author Jian
 *
 */
public class MixDBSearcher extends SimpleDBSearcher{
	public int minRank; //how deep we go down the initial list to get the first spectrum
	public int maxRank; //how deep we go down the initial list to get the second spectrum
	public String mixScorerPath="/resources/yeast_simmix_alpha_generic_8_25.o";
	public SpectrumComparator mixScorer;
	public String precursorsFile;
	public String ptmFile;
	public double windowWidth=25; //the isolation windowWidth
	public Map<Integer, double[]> precursorList;
	public List<PTM> ptms;
	public int maxPTM=2;
	public List<PTM[]> ptmList;
	public int minScan=0;
	public int maxScan=Integer.MAX_VALUE;
	
	public MixDBSearcher(String dbPath, String spectrumFile, double tolerance, String outputFile){
		this(dbPath, spectrumFile, tolerance, 0.5, outputFile);
	}
	
	public MixDBSearcher(String dbPath, String spectrumFile, double tolerance, double fragmentTolerance, String outputFile){
		this(dbPath, spectrumFile, tolerance, fragmentTolerance, outputFile, 0, Integer.MAX_VALUE);
	}
		
	public MixDBSearcher(String dbPath, String spectrumFile, double tolerance, double fragmentTolerance, String outputFile, int minScan, int maxScan){
		super(dbPath, spectrumFile);
		this.parentTolerance = tolerance;
		this.outputFile = outputFile;
		this.fragmentTolerance = fragmentTolerance;
		this.minScan = minScan;
		this.maxScan = maxScan;
	}
	
	
	/**
	 * Search query spectra against database, look for the best pair of peptie matches to each spectrum
	 */
	public void searchDB(){
		initialize();
		//MZXMLReader reader = new MZXMLReader(this.spectrumFile);
		//MZXMLReader reader = new SortedMZXMLReader(this.spectrumFile);
		Iterator reader = null;
		if(this.spectrumFile.endsWith(".mgf")){
			reader = new SortedSpectrumReader(this.spectrumFile);
		}
		
		if(this.spectrumFile.toLowerCase().endsWith("mzxml")){
			reader = new SortedMZXMLReader(this.spectrumFile, this.minScan, this.maxScan);
		}
		
		RankBaseScoreLearner pcomp = RankBaseScoreLearner.loadComparatorLocal(this.singleScorer);
		SimpleProbabilisticScorer scorer3 = new SimpleProbabilisticScorer(pcomp);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int counter = 0;
		System.out.println("start searching");
		BufferedWriter out = initOutput(this.outputFile);
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			if(s.getPeak().size() < 10){
				continue;
			}
			//s.windowFilterPeaks(topPeakKept, windowWidth);
			s.windowFilterPeaks(12, 25);
			s.computePeakRank();
//			String[] peptides = s.peptide.split(" & ");
//			Peptide p1 = new Peptide(peptides[0]);
//			Peptide p2 = new Peptide(peptides[1]);
//			System.out.println(s.scanNumber + "\t" + s.spectrumName + "\t" + p1 + "\t" + p2 + "\t" + p1.getParentmass() + "\t" + p2.getParentmass());
			ArraySpectrum a = ArraySpectrum.getRankSpectrum(s);
			List<Spectrum> cands = new ArrayList();
			for(int c = minCharge; c <= maxCharge; c++){
				double pm = s.parentMass*c-Mass.WATER-c*Mass.PROTON_MASS;
				double tolerance = this.parentTolerance*c;
				cands.addAll(this.theoDB.getCandidates(s, this.parentTolerance));
			}
			//System.out.println("Number of candidates: " + cands.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cands, this.comp, this.mixScorer);
			searcher.spectrumFile = this.spectrumFile;
			searcher.bw = out;
			searcher.setSingleScorer(this.comp);
			searcher.bestArrayCandidates(a, 1);
			counter++;
		}
		try{
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		System.out.println("matching " + counter + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	/**
	 * Same as searchDB(), but extends it to look for Mix number of peptides per spectrum
	 * @param Mix
	 */
	public void searchDB(int Mix){
		initialize();
		//MZXMLReader reader = new MZXMLReader(this.spectrumFile);
		//MZXMLReader reader = new SortedMZXMLReader(this.spectrumFile);
		Iterator reader = null;
		if(this.spectrumFile.endsWith(".mgf")){
			reader = new SortedSpectrumReader(this.spectrumFile);
		}
		
		if(this.spectrumFile.toLowerCase().endsWith("mzxml")){
			reader = new SortedMZXMLReader(this.spectrumFile, this.minScan, this.maxScan);
		}
		
		RankBaseScoreLearner pcomp = RankBaseScoreLearner.loadComparatorLocal(this.singleScorer);
		SimpleProbabilisticScorer scorer3 = new SimpleProbabilisticScorer(pcomp);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int counter = 0;
		System.out.println("start searching");
		BufferedWriter out = initOutput(this.outputFile);
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			if(s.getPeak().size() < 10){
				continue;
			}
			if(s.scanNumber != 27){
				//continue;
			}
			//s.windowFilterPeaks(topPeakKept, windowWidth);
			s.windowFilterPeaks(10, 25);
			s.computePeakRank();
//			String[] peptides = s.peptide.split(" & ");
//			Peptide p1 = new Peptide(peptides[0]);
//			Peptide p2 = new Peptide(peptides[1]);
//			System.out.println(s.scanNumber + "\t" + s.spectrumName + "\t" + p1 + "\t" + p2 + "\t" + p1.getParentmass() + "\t" + p2.getParentmass());
			ArraySpectrum a = ArraySpectrum.getRankSpectrum(s);
			List<Spectrum> cands = new ArrayList();
			minCharge = 1; maxCharge = 1;
			for(int c = minCharge; c <= maxCharge; c++){
				double pm = s.parentMass*c-Mass.WATER-c*Mass.PROTON_MASS;
				double tolerance = this.parentTolerance*c;
				cands.addAll(this.theoDB.getCandidates(s, this.parentTolerance));
			}
			System.out.println(s.spectrumName + "\t" +  s.parentMass + "\t" + s.charge + "\tNumber of candidates: " + cands.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cands, this.comp, this.mixScorer);
			searcher.spectrumFile = this.spectrumFile;
			searcher.bw = out;
			searcher.setSingleScorer(this.comp);
			searcher.checkMixtureRanks(a);
			//searcher.bestArrayCandidates(a, 1, Mix);
			counter++;
		}
		try{
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		System.out.println("matching " + counter + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	/**
	 * This search allow one to specify a list of precursors to search for each query spectrum
	 */
	public void searchWithMultiPrecursors(){
		initialize();
		Iterator reader = null;
		if(this.spectrumFile.endsWith(".mgf")){
			reader = new SortedSpectrumReader(this.spectrumFile);
		}
		
		if(this.spectrumFile.toLowerCase().endsWith("mzxml")){
			reader = new SortedMZXMLReader(this.spectrumFile, this.minScan, this.maxScan);
		}
		
		RankBaseScoreLearner pcomp = RankBaseScoreLearner.loadComparatorLocal(this.singleScorer);
		SimpleProbabilisticScorer scorer3 = new SimpleProbabilisticScorer(pcomp);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int counter = 0;
		System.out.println("start searching");
		System.out.println("Using mix_model " + this.mixScorerPath);
		BufferedWriter out = initOutput(this.outputFile);
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			if(s.getPeak().size() < 10 || s.scanNumber < this.minScan || s.scanNumber >= this.maxScan){
				continue;
			}
			//s.windowFilterPeaks(topPeakKept, windowWidth);
			s.windowFilterPeaks(12, 25);
			s.computePeakRank();
//			String[] peptides = s.peptide.split(" & ");
//			Peptide p1 = new Peptide(peptides[0]);
//			Peptide p2 = new Peptide(peptides[1]);
//			System.out.println(s.scanNumber + "\t" + s.spectrumName + "\t" + p1 + "\t" + p2 + "\t" + p1.getParentmass() + "\t" + p2.getParentmass());
			ArraySpectrum a = ArraySpectrum.getRankSpectrum(s);
			List<Spectrum> cands = new ArrayList();
			
			int scan = s.scanNumber;
			double[] precursors = null;
			if(this.precursorList.get(scan) == null){
				continue;
			}else{
				precursors = this.precursorList.get(scan);
			}
			
			System.out.println("Number of detected precursors: " + (precursors.length/2));
			for(int i = 0; i < precursors.length-1; i=i+2){
				//s.parentMass = precursors[i];
				//s.charge = (int)precursors[i+1];
				//System.out.println("precursor: " + s.parentMass + "\t" + s.charge);
				int charge = (int)precursors[i+1];
				cands.addAll(this.theoDB.getCandidates(s, this.windowWidth, precursors[i], this.parentTolerance, charge));
				
				if(this.parentTolerance < 1.0){
					double c13Offset = Mass.C13 - Mass.C12;
					cands.addAll(this.theoDB.getCandidates(s, this.windowWidth, precursors[i]+c13Offset/charge, this.parentTolerance, charge));
					cands.addAll(this.theoDB.getCandidates(s, this.windowWidth, precursors[i]-c13Offset/charge, this.parentTolerance, charge));
					//cands.addAll(this.theoDB.getCandidates(s, this.windowWidth, precursors[i]+c13Offset/3, this.parentTolerance, (int)precursors[i+1]));
					//cands.addAll(this.theoDB.getCandidates(s, this.windowWidth, precursors[i]-c13Offset/3, this.parentTolerance, (int)precursors[i+1]));
				}
			}
			System.out.println("Number of candidates: " + cands.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cands, this.comp, this.mixScorer);
			searcher.spectrumFile = this.spectrumFile;
			searcher.bw = out;
			searcher.setSingleScorer(this.comp);
			searcher.bestArrayCandidates(a, 1);
			counter++;
		}
		try{
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		System.out.println("matching " + counter + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	/**
	 * This search add PTM support for search.  PTMs are input as a parameter files
	 */
	public void searchDBWithPTM(){
		initialize();
		Iterator reader = null;
		if(this.spectrumFile.endsWith(".mgf")){
			reader = new SortedSpectrumReader(this.spectrumFile);
		}
		
		if(this.spectrumFile.toLowerCase().endsWith("mzxml")){
			reader = new SortedMZXMLReader(this.spectrumFile);
			//reader = new MZXMLReader(this.spectrumFile)
		}
		
		RankBaseScoreLearner pcomp = RankBaseScoreLearner.loadComparatorLocal(this.singleScorer);
		SimpleProbabilisticScorer scorer3 = new SimpleProbabilisticScorer(pcomp);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int counter = 0;
		System.out.println("start searching");
		BufferedWriter out = initOutput(this.outputFile);
		while(reader.hasNext()){
			Spectrum s = (Spectrum)reader.next();
			if(s.getPeak().size() < 10){
				continue;
			}
			if(s.scanNumber != 9552){
				//continue;
			}
			//s.windowFilterPeaks(topPeakKept, windowWidth);
			s.windowFilterPeaks(10, 25);
			s.computePeakRank();
			ArraySpectrum a = ArraySpectrum.getRankSpectrum(s);
			List<Spectrum> cands = new ArrayList();
			for(int c = minCharge; c <= maxCharge; c++){
				double pm = s.parentMass*c-Mass.WATER-c*Mass.PROTON_MASS;
				double tolerance = this.parentTolerance*c;
				cands.addAll(this.theoDB.getCandidates(s, this.parentTolerance, this.ptmList));
			}
			System.out.println("Number of candidates: " + cands.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cands, this.comp, this.mixScorer);
			searcher.spectrumFile = this.spectrumFile;
			searcher.bw = out;
			searcher.setSingleScorer(this.comp);
			searcher.bestArrayCandidates(a, 1);
			counter++;
		}
		try{
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		System.out.println("matching " + counter + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public BufferedWriter initOutput(String outFile){
		BufferedWriter bw=null;
		if(this.outputFile != null){
			try{
				bw = new BufferedWriter(new FileWriter(this.outputFile));
			}catch(IOException ioe){
				System.err.println(ioe.getMessage());
				ioe.printStackTrace();
			}
		}else{
			bw = new BufferedWriter(new OutputStreamWriter(System.out));
		}
		return bw;
	}
	
	/**
	 * If spectrum has multiple precursors and such information can be input
	 * to MixDB using a file that has a list of precursors for each spectrum
	 * @param precursorsFile
	 */
	public Map<Integer, double[]> parsePrecursorList(String precursorsFile){
		List<String> lines = Utils.FileIOUtils.createListFromFile(precursorsFile);
		Map<Integer, List<Double>> precursorList = new HashMap<Integer, List<Double>>();
		for(int i = 0; i < lines.size(); i++){
			String[] tokens = lines.get(i).split("\\t");
			int scan = Integer.parseInt(tokens[0]);
			List<Double> precursors;
			if(precursorList.containsKey(scan)){
				precursors = precursorList.get(scan);
			}else{
				precursors = new ArrayList<Double>();
				precursorList.put(scan, precursors);
			}
			precursors.add(Double.parseDouble(tokens[1]));
			precursors.add(Double.parseDouble(tokens[2]));
		}
		Map<Integer, double[]> precursorList2 = new HashMap<Integer, double[]>();
		for(Iterator<Integer> it = precursorList.keySet().iterator(); it.hasNext();){
			int scan = it.next();
			List<Double> precursors = precursorList.get(scan);
			double[] precursors2 = new double[precursors.size()];
			for(int i = 0; i < precursors.size(); i++){
				precursors2[i] = precursors.get(i);
			}
			precursorList2.put(scan, precursors2);
		}
		System.out.println("Finish getting precursors info, precursor map size: " + precursorList2.keySet().size());
		return precursorList2;
	}
	
	/**
	 * Parse PTM files to generate the possible ptms allowed for this search
	 * @param precursorsFile
	 * @return
	 */
	public List<PTM[]> parsePTMs(String ptmFile){
		List<PTM[]> ptmList = new ArrayList();
		List<String> lines = Utils.FileIOUtils.createListFromFile(ptmFile);
		int maxPtm = 0;
		this.ptmFile = ptmFile;
		this.ptms = new ArrayList<PTM>();
//		this.ptms.add(new PTM(28.03, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(32.06, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(4.01, new int[]{PTM.CTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(16.01, new int[]{PTM.ANYPOSITION}, new char[]{'M'}));
		for(int i = 0; i < lines.size(); i++){
			String[] line = lines.get(i).split(",");
			if(line[0].startsWith("#")){
				continue;
			}else if(line[0].equals("maxPTM")){
				maxPtm = Integer.parseInt(line[1]);
			}else{
				PTM ptm = this.parsePTM(line);
				this.ptms.add(ptm);
			}
		}
		this.ptmList = PTM.generatePTMList(this.ptms, maxPtm);	
		ptmList.add(new PTM[]{});
		ptmList.addAll(this.ptmList);
		this.ptmList = ptmList;
		System.out.println("Finish getting ptms info, ptm size: " + this.ptmList.size());
		return ptmList;
	}
	
	private PTM parsePTM(String[] tokens){
		double PtmMass = Double.parseDouble(tokens[0]);
		String residues = tokens[1];
		char[] targetedRes = new char[residues.length()];
		int[] targetedPos = new int[residues.length()];
		int pos=0;
		for(int i = 0; i < residues.length(); i++){
			if(residues.charAt(i) == '*'){
				targetedRes[i] = PTM.ANYRESIDUE;
			}else{
				targetedRes[i] = residues.charAt(i);
			}
		}
		
		if(tokens[2].equals("N-term")){
			pos = PTM.NTERM;
		}else if(tokens[2].equals("C-term")){
			pos = PTM.CTERM;
		}else if(tokens[2].equals("any")){
			pos = PTM.ANYPOSITION;
		}else{
			pos = Integer.parseInt(tokens[2]);
		}
		targetedPos = new int[]{pos};
		return new PTM(PtmMass, targetedPos, targetedRes);
	}
	
	public void initialize(){
		this.db = new DatabaseIndexer(this.dbPath);
		this.theoDB = new TheoSpectrumGenerator(this.db);
		RankBaseScoreLearner pComp = RankBaseScoreLearner.loadComparatorLocal(this.singleScorer);
		this.comp = new SimpleProbabilisticScorer(pComp);
		((SimpleProbabilisticScorer)this.comp).setMinMatchedPeak(0);
		this.comp = ArraySpectrumComparator.loadStandardComparator(singleScorer);
		this.mixScorer = ArraySpectrumComparator.loadMixtureComparator(mixScorerPath);
		((ArraySpectrumComparator)this.comp).massTolerance = this.fragmentTolerance;
		((ArraySpectrumComparator)this.mixScorer).massTolerance = this.fragmentTolerance;
		if(this.precursorsFile != null){
			this.precursorList = this.parsePrecursorList(this.precursorsFile); //if the optinal precursor list is inputed we also parse this information 
		}
		if(this.ptmFile != null){
			this.parsePTMs(this.ptmFile);
		}
	}
	
	@Override
	public TreeSet topNCandidates(Spectrum s, Collection db, int topN) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public static void initSearch(String[] args){
		Parser parser = new GnuParser();
    	Options options = createOptions();
    	try{
    		CommandLine line = parser.parse( options, args);
    		if(line.hasOption("help")){
    			HelpFormatter formatter = new HelpFormatter();
    			formatter.printHelp( "MixDB", options );
    			return;
    	     }	
    		
    	    String spectFile = line.getOptionValue("s");
 	        String dbFile = line.getOptionValue("d");
 	        String outFile = line.getOptionValue("o");
 	        double pmTol = Double.parseDouble(line.getOptionValue("p"));
 	        double fmTol = Double.parseDouble(line.getOptionValue("f"));
 	        if(fmTol > 0.5){
 	        	System.err.println("Warning : fragment mass tolerance should be less than 0.5, using 0.5 as fragment mass tolerance");
 	        	fmTol = 0.5;
 	        }
 	        MixDBSearcher searcher = new MixDBSearcher(dbFile, spectFile, pmTol, fmTol, outFile);
    	    if(line.hasOption("m")){
    	    	searcher.mixScorerPath = line.getOptionValue("m");
    	    }
 	        if(line.hasOption("minScan")){
 	        	searcher.minScan = Integer.parseInt(line.getOptionValue("minScan"));
 	        }
 	        if(line.hasOption("maxScan")){
 	        	searcher.maxScan = Integer.parseInt(line.getOptionValue("maxScan"));
 	        }
    	    if(line.hasOption("precursorList")){
    	    	searcher.precursorsFile	= line.getOptionValue("precursorList");
    	    	searcher.searchWithMultiPrecursors();
    	    }else if(line.hasOption("ptmList")){
    	    	searcher.ptmFile = line.getOptionValue("ptmList");
    	    	searcher.searchDBWithPTM();
    	    }else{ 	   	    
    	    	searcher.searchDB();
    	    }
    		
    	}catch( ParseException exp ) {
	        // oops, something went wrong
	        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
	        HelpFormatter formatter = new HelpFormatter();
        	formatter.printHelp( "MixDB", options);
	    }catch( NumberFormatException e){
	    	System.err.println("Mass tolerance or scan parameters must be a valid number");
	    	e.getMessage();
	    	e.printStackTrace();
		}catch(Exception e){
	    	System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(-1);
	    }
    	
   }
	
	private static Options createOptions(){
		Option help = new Option("help", "help");
		Option fastaFile = OptionBuilder.withArgName( "fasta file" )
                .hasArg()
                .withDescription(  "sequence database file in fasta" )
                .isRequired()
                .create( "d" );
		
		Option spectrumFile = OptionBuilder.withArgName( "spectrum file" )
                .hasArg()
                .withDescription(  "spectrum file" )
                .isRequired()
                .create( "s" );
		
		Option pmTol = OptionBuilder.withArgName( "parent mass tolerance" )
                .hasArg()
                .withDescription(  "A window of +/- K Da around observed mass to consider peptide candidate matches" )
                .isRequired()
                .create( "p" );
		
		Option fmTol = OptionBuilder.withArgName( "fragment mass tolerance" )
                .hasArg()
                .withDescription(  "A window of +/- K Da around observed mass to consider fragment ions matches" )
                .isRequired()
                .create( "f" );
		
		
		Option minScan = OptionBuilder.withArgName( "min Scan" )
                .hasArg()
                .withDescription(  "minimum scan to search" )
                .create( "minScan" );
		
		Option maxScan = OptionBuilder.withArgName( "max Scan" )
                .hasArg()
                .withDescription(  "maximum scan to search" )
                .create( "maxScan" );
		
		Option outFile = OptionBuilder.withArgName( "output file" )
                .hasArg()
                .withDescription(  "output for search result" )
                .isRequired()
                .create( "o" );
		
		Option precursorList = OptionBuilder.withArgName( "precursor list" )
                .hasArg()
                .withDescription(  "specificy list of precursors for each spectrum" )
                .create( "precursorList" );
		
		Option ptmList = OptionBuilder.withArgName( "ptm list" )
                .hasArg()
                .withDescription(  "specificy list of ptm to search" )
                .create( "precursorList" );
		
		Option scoreModel = OptionBuilder.withArgName( "model" )
                .hasArg()
                .withDescription(  "scoring model file" )
                .create( "m" );
		
		Options options = new Options();
		options.addOption(help);
		options.addOption(outFile);
		options.addOption(fastaFile);
		options.addOption(spectrumFile);
		options.addOption(pmTol);
		options.addOption(fmTol);
		options.addOption(minScan);
		options.addOption(maxScan);
		options.addOption(scoreModel);
		return options;
	}
	
	public static void main(String[] args){
		//args[0] = "../mixture_linked/database/Ecoli_genome_plusDecoy.fasta";
		//args[1] = "../../Downloads/ACG_Nuno-selected/New_Swath_3mu/SWATH_3amu_bin16_01.mzXML";
		//args[2] = "3.0";
		//args[3] = "0.05";
		//args[4] = "../mixture_linked/ACG_14344_hardklorPrecursorsListmin07.txt";
		//args[5] = "../mixture_linked/Mod_mixdb.txt";
		//args[5] ="../mixture_linked/testmixdb.txt";
		//args[6] = "1";
		//args[7] ="750000";
		//args[8] = "mixtures_TOF_alpha01-10_models.o";
		initSearch(args);
	}
}
