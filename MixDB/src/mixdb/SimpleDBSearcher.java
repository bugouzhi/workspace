package mixdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import org.Spectrums.LazyEvaluatedSpectrum;
import org.Spectrums.MZXMLReader;
import org.Spectrums.Mass;
import org.Spectrums.Peptide;
import org.Spectrums.PeptideLite;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.SimpleProbabilisticScorer;
import org.Spectrums.SortedMZXMLReader;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumComparator;
import org.Spectrums.SpectrumLibSearcher;
import org.Spectrums.TheoreticalSpectrum;

import sequences.FastaSequence;

/**
 * An implementation of the DBSearcher interface
 * that aim at searching single-peptide per query spectrum
 * (as oppotion to MixDBSearcher which search multi-peptide matches to a query spectrum).
 * @author Jian Wang
 *
 */
public class SimpleDBSearcher implements DBSearcher{
	public DatabaseIndexer db;
	public TheoSpectrumGenerator theoDB;
	public SpectrumComparator comp;
	public TreeSet topCandidates;
	public String spectrumFile;
	public String dbPath;
	public String singleScorer="/resources/yeast_single_model_realannotated_win10_25.o";
	public double parentTolerance = 1.5;
	public double fragmentTolerance = 0.05;
	public int topPeakKept = 10;
	public int windowWidth = 25;
	public int minCharge = 2;
	public int maxCharge = 3;
	public boolean matchCharge=false;
	public String outputFile="./mixdb.out"; //default output

	public SimpleDBSearcher(String dbPath, String spectrumFile){
		this.spectrumFile = spectrumFile;
		this.dbPath = dbPath;
	}
	
	
	public void searchDB(){
		initialize();
		//MZXMLReader reader = new MZXMLReader(this.spectrumFile);
		//MZXMLReader reader = new SortedMZXMLReader(this.spectrumFile);
		Iterator<Spectrum> reader = new SortedSpectrumReader(this.spectrumFile);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int counter=0;
		while(reader.hasNext()){
			Spectrum s = reader.next();
			if(s.getPeak().size() < 10 || s.scanNumber != 26694){
				//continue;
			}
			s.windowFilterPeaks(topPeakKept, windowWidth);	
			s.computePeakRank();
			ArraySpectrum a = ArraySpectrum.getRankSpectrum(s);
			List<Spectrum> cands = new ArrayList();
			for(int c = minCharge; c <= maxCharge; c++){
				double pm = s.parentMass*c-Mass.WATER-c*Mass.PROTON_MASS;
				double tolerance = this.parentTolerance*c;
				cands.addAll(this.theoDB.getCandidates(s, this.parentTolerance));
			}
			System.out.println("Scan\t" + s.scanNumber + "\tNumber of candidates:\t" + cands.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cands, this.comp);
			searcher.spectrumFile = this.spectrumFile;
			searcher.bestSpectra(a, 2);
			counter++;
		}
		System.out.println("matching " + counter + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000);
	}
	
	public List getTheoreticalSpectrum(List<PeptideLite> peps, Spectrum s, int charge){
		List candidates = new ArrayList();
		FastaSequence seq = this.db.getSeq();
		for(int i = 0; i < peps.size(); i++){
			PeptideLite pep = peps.get(i);
			Peptide p = new Peptide(seq.getSubsequence(pep.getBeginInd(), pep.getEndInd()+1)+"."+charge);
			//TheoreticalSpectrum th = new TheoreticalSpectrum(p, new String[]{"b"}, new String[]{"y"});
			LazyEvaluatedSpectrum th = new LazyEvaluatedSpectrum(p);
			candidates.add(th);
		}
		return candidates;
	}
	
	
	public List getArrayTheoreticalSpectrum(List<PeptideLite> peps, Spectrum s, int charge){
		List candidates = new ArrayList();
		FastaSequence seq = this.db.getSeq();
		for(int i = 0; i < peps.size(); i++){
			PeptideLite pep = peps.get(i);
			Peptide p = new Peptide(seq.getSubsequence(pep.getBeginInd(), pep.getEndInd()+1)+"."+charge);
			ArrayTheoreticalSpectrum th = (ArrayTheoreticalSpectrum)TheoreticalSpectrumFactory.getTheoSpectrumX(p.getPeptide(), p.getCharge(), 
					TheoreticalSpectrumFactory.standardTypeMap, 
					TheoreticalSpectrumFactory.standardIonMap);
			
			candidates.add(th);
		}
		return candidates;
	}
	
	public void initialize(){
		this.db = new DatabaseIndexer(this.dbPath);
		this.theoDB = new TheoSpectrumGenerator(this.db);
		RankBaseScoreLearner pComp = RankBaseScoreLearner.loadComparatorLocal(this.singleScorer);
		this.comp = new SimpleProbabilisticScorer(pComp);
		((SimpleProbabilisticScorer)this.comp).setMinMatchedPeak(0);
		this.comp = ArraySpectrumComparator.loadStandardComparator(singleScorer);
		((ArraySpectrumComparator)this.comp).massTolerance = this.fragmentTolerance;
	}
	
	@Override
	public TreeSet topNCandidates(Spectrum s, Collection db, int topN) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public static void main(String[] args){
		args[0] = "../mixture_linked/database/Human_allproteins.fasta";
		args[1] = "../mixture_linked/msdata/Training/MSGFDB_Tryp_7.mgf";
		SimpleDBSearcher searcher = new SimpleDBSearcher(args[0], args[1]);
		searcher.searchDB();
	}
	
	
	
}
