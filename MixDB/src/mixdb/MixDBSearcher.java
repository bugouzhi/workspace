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
import Utils.FileIOUtils;

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
	public int minScan;
	public int maxScan;
	
	public MixDBSearcher(String dbPath, String spectrumFile, double tolerance, String outputFile){
		super(dbPath, spectrumFile);
		this.parentTolerance = tolerance;
		this.outputFile = outputFile;
		this.fragmentTolerance = 0.5;
	}
	
	public MixDBSearcher(String dbPath, String spectrumFile, double tolerance, double fragmentTolerance, String outputFile){
		super(dbPath, spectrumFile);
		this.parentTolerance = tolerance;
		this.outputFile = outputFile;
		this.fragmentTolerance = fragmentTolerance;
	}
	
	
	public void searchDB(){
		initialize();
		//MZXMLReader reader = new MZXMLReader(this.spectrumFile);
		//MZXMLReader reader = new SortedMZXMLReader(this.spectrumFile);
		Iterator reader = null;
		if(this.spectrumFile.endsWith(".mgf")){
			reader = new SortedSpectrumReader(this.spectrumFile);
		}
		
		if(this.spectrumFile.toLowerCase().endsWith("mzxml")){
			reader = new SortedMZXMLReader(this.spectrumFile);
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
			for(int c = minCharge; c <= maxCharge; c++){
				double pm = s.parentMass*c-Mass.WATER-c*Mass.PROTON_MASS;
				double tolerance = this.parentTolerance*c;
				cands.addAll(this.theoDB.getCandidates(s, this.parentTolerance));
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
	
	
	public void searchWithMultiPrecursors(){
		initialize();
		Iterator reader = null;
		if(this.spectrumFile.endsWith(".mgf")){
			reader = new SortedSpectrumReader(this.spectrumFile);
		}
		
		if(this.spectrumFile.toLowerCase().endsWith("mzxml")){
			reader = new SortedMZXMLReader(this.spectrumFile);
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
			
			int scan = s.scanNumber;
			double[] precursors = null;
			if(this.precursorList.get(scan) == null){
				continue;
			}else{
				precursors = this.precursorList.get(scan);
			}
			
			for(int i = 0; i < precursors.length-1; i=i+2){
				//s.parentMass = precursors[i];
				//s.charge = (int)precursors[i+1];
				//System.out.println("precursor: " + s.parentMass + "\t" + s.charge);
				cands.addAll(this.theoDB.getCandidates(s, this.windowWidth, precursors[i], this.parentTolerance, (int)precursors[i+1]));
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
	
	public void searchDBWithPTM(){
		initialize();
		Iterator reader = null;
		if(this.spectrumFile.endsWith(".mgf")){
			reader = new SortedSpectrumReader(this.spectrumFile);
		}
		
		if(this.spectrumFile.toLowerCase().endsWith("mzxml")){
			reader = new SortedMZXMLReader(this.spectrumFile);
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
		this.ptmList = PTM.generatePTMList(this.ptms, 3);	
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
	
	public static void main(String[] args){
		args[0] = "../mixture_linked/database/RpalTIE1.fasta";
		args[1] = "../mixture_linked/TIE1_IPTL_sample_mixtures.mgf";
		args[2] = "0.05";
		args[3] = "0.5";
		args[4] = "../mixture_linked/ACG_14348_hardklorPrecursorsList.txt";
		args[5] = "../mixture_linked/Mod_mixdb.txt";
		args[6] ="../mixture_linked/testmixdbPtm.txt";
		if(args.length != 6 && args.length != 7){
			System.out.println("usage: java -Xmx2000M -jar MixDB.jar <database> <spectraFile> <parentmass tolerance> <fragment mass tolerance> <modification file> <outfile>");
			System.out.println("   or: java -Xmx2000M -jar MixDB.jar <database> <spectraFile> <parentmass tolerance> <fragment mass tolerance> <precursor list> <modification file> <outfile>");
			return;
		}
		
		MixDBSearcher searcher;
		if(args.length == 7 && args[4].length() > 1){
			searcher = new MixDBSearcher(args[0], args[1], 
					Double.parseDouble(args[2]), Double.parseDouble(args[3]), args[6]);
			searcher.precursorsFile = args[4];
		}else{
			searcher = new MixDBSearcher(args[0], args[1], 
				Double.parseDouble(args[2]), Double.parseDouble(args[3]), args[4]);
		}
		searcher.spectrumFile = args[1];
		searcher.ptmFile = args[5];
		System.out.println("Start running");
		try{
			if(args.length == 6){
				searcher.searchWithMultiPrecursors();
			}else if(args.length == 7){
				searcher.searchDBWithPTM();
			}else{
				searcher.searchDB();
			}
		}catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
