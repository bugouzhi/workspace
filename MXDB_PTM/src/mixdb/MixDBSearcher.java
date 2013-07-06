package mixdb;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import org.Spectrums.LazyEvaluatedSpectrum;
import org.Spectrums.MZXMLReader;
import org.Spectrums.Mass;
import org.Spectrums.MixturePeakScoreLearner;
import org.Spectrums.Peptide;
import org.Spectrums.PeptideLite;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.SimpleProbabilisticScorer;
import org.Spectrums.SortedMZXMLReader;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumComparator;
import org.Spectrums.SpectrumLibSearcher;

import SeqDB.DatabaseIndexer;

public class MixDBSearcher extends SimpleDBSearcher{
	public int minRank; //how deep we go down the initial list to get the first spectrum
	public int maxRank; //how deep we go down the initial list to get the second spectrum
	public String mixScorerPath="/resources/yeast_simmix_alpha_generic_8_25.o";
	public SpectrumComparator mixScorer;
	
	public MixDBSearcher(String dbPath, String spectrumFile, double tolerance, String outputFile){
		super(dbPath, spectrumFile);
		this.parentTolerance = tolerance;
		this.outputFile = outputFile;
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
			s.windowFilterPeaks(8, 25);
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
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cands, this.comp, this.mixScorer);
			searcher.bw = out;
			searcher.setSingleScorer(this.comp);
			searcher.bestArrayCandidates(a, 10);
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
	
	public void initialize(){
		this.db = new DatabaseIndexer(this.dbPath);
		this.theoDB = new TheoSpectrumGenerator(this.db);
		RankBaseScoreLearner pComp = RankBaseScoreLearner.loadComparatorLocal(this.singleScorer);
		this.comp = new SimpleProbabilisticScorer(pComp);
		((SimpleProbabilisticScorer)this.comp).setMinMatchedPeak(0);
		this.comp = ArraySpectrumComparator.loadStandardComparator(singleScorer);
		this.mixScorer = ArraySpectrumComparator.loadMixtureComparator(mixScorerPath);
		
	}
	
	@Override
	public TreeSet topNCandidates(Spectrum s, Collection db, int topN) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public static void main(String[] args){
		args[0] = "../mixture_linked/database/Human_allproteins_plusDecoy.fasta";
		args[1] = "../mixture_linked/human_heck_simmix.mgf";
		args[2] = "3.0";
		args[3] ="../mixture_linked/human_simmix_mixdb.txt";
		MixDBSearcher searcher = new MixDBSearcher(args[0], args[1], Double.parseDouble(args[2]), args[3]);
		searcher.searchDB();
	}
}
