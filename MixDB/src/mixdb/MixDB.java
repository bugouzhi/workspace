package mixdb;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.List;

import org.Spectrums.LazyEvaluatedSpectrum;
import org.Spectrums.MZXMLReader;
import org.Spectrums.Mass;
import org.Spectrums.MixturePeakScoreLearner;
import org.Spectrums.MixtureSpectrumScorer;
import org.Spectrums.PeakComparator;
import org.Spectrums.Peptide;
import org.Spectrums.PeptideLite;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.SimplePeakComparator;
import org.Spectrums.SimpleProbabilisticScorer;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumComparator;
import org.Spectrums.SpectrumLib;
import org.Spectrums.SpectrumLibSearcher;
import org.Spectrums.TheoreticalSpectrum;

import sequences.FastaSequence;

public class MixDB {
	private String spectrumFile;
	private String dbPath;
	private String mixScorer="../mixture_linked/yeast_simmix_alpha_generic_12_25.o";
	private String singleScorer="../mixture_linked/yeast_single_model_realannotated_win10_25.o";
	private double parentTolerance = 3.0;
	private double fragmentTolerance = 0.5;
	private int topPeakKept = 8;
	private int windowWidth = 25;
	private int minCharge = 2;
	private int maxCharge = 3;
	private boolean matchCharge=false;
	private DatabaseIndexer db;
	
	public MixDB(String dbPath, String spectrumFile){
		this.dbPath = dbPath;
		this.spectrumFile = spectrumFile;
		this.db = new DatabaseIndexer(this.dbPath);
		this.db.indexDatabase();
	}
	
	public void runMixDB(){
		String probFile = "..//mixture_linked//data//IonsScore.txt";
		String noiseModel = "..//mixture_linked//data//NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer filter = new SimpleProbabilisticScorer(peakscorer);
		filter.setIncludeNoise(false);
		RankBaseScoreLearner peakscorer2 = RankBaseScoreLearner.loadComparator(singleScorer);
		SimpleProbabilisticScorer scorer1 = new SimpleProbabilisticScorer(peakscorer2);
		scorer1.setMinMatchedPeak(0);
		PeakComparator peakscorer3 = MixturePeakScoreLearner.loadComparator(mixScorer);
		MixtureSpectrumScorer scorer = new MixtureSpectrumScorer(peakscorer3);
		MZXMLReader iterator = new MZXMLReader(spectrumFile);
		SpectrumComparator scorer3 = ArraySpectrumComparator.loadStandardComparator(singleScorer);
		long start = (new GregorianCalendar()).getTimeInMillis();
		int counter = 1;
		for(; iterator.hasNext();){
			Spectrum s = iterator.next();
			if(s.getPeak().size() < 10){
				continue;
			}
			s.windowFilterPeaks(topPeakKept, windowWidth);
			s.computePeakRank();
			ArraySpectrum a = ArraySpectrum.getRankSpectrum(s);
			if(s.scanNumber < 1){
				s.scanNumber = counter;
			}
			if(s.scanNumber != 6660){
				//continue;
			}
			List<Spectrum> cands = new ArrayList();
			for(int c = minCharge; c <= maxCharge; c++){
				double pm = s.parentMass*c-Mass.WATER-c*Mass.PROTON_MASS;
				double tolerance = this.parentTolerance*c;
				//System.out.println("looking for mass: " + pm + "tolerance: " + tolerance);
				List<PeptideLite> pep = this.db.getPeptides(pm-tolerance, 
						pm+tolerance, tolerance);
				//cands.addAll(getTheoreticalSpectrum(pep, s, c));
				cands.addAll(getArrayTheoreticalSpectrum(pep, s, c));
			}
			//System.out.println("Scan: " + s.scanNumber + "\tnumber of peps: " + cands.size());
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(cands, scorer3, scorer);
			searcher.setSingleScorer(scorer3);
			//searcher.bestCandidates(s, 10);
			searcher.bestCandidates(a, 10);
			counter++;
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
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
			ArrayTheoreticalSpectrum th = TheoreticalSpectrumFactory.getArrayTheoSpectrum(p, 
					TheoreticalSpectrumFactory.standardTypeMap, 
					TheoreticalSpectrumFactory.standardIonMap);
			candidates.add(th);
		}
		return candidates;
	}

	
	public static void main(String[] args){
		args[0] = "../mixture_linked/database/yeast_proteins_plusDecoy.fasta";
		args[1] = "../mixture_linked/mixtures100000_alpha1.0.mgf";
		MixDB mixdb = new MixDB(args[0], args[1]);
		mixdb.runMixDB();
	}
	
}
