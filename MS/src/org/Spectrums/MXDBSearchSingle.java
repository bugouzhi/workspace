package org.Spectrums;
import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import mixgf.PRMSpectrum;
import mixgf.PRMSpectrumComparator;

public class MXDBSearchSingle {
	public double windowWidth = 25;
	public int topPeaksKept = 8;
	//candidate filtering
	public double parentMassTolerance = 2;
	public int minMatchedPeak = 3; //for filtering
	public int minContinuousMatch = 1; // for filtering
	public int topFilteringPeak = 20; //for filtering
	public int topCandidateKept = 100;	
	public double fragmentMassTolerance = 0.5;
	//input and training file
	public String queryFile;
	public String peptidesFile;
	private String training = "..\\MSPLib\\Lib\\ecoli.msp";
	private String mixtureTraining = "";
	
	//peptide table
	TheoreticalCandidateSpectrumFactory factory;
	Iterator<Spectrum> queryIterator;
	SpectrumComparator comp;
	PeakComparator pComp;
	TreeMap<Double, Spectrum> topCandidates;
	private double currentMax;
	
	public static void  buildPeptideMap(String peptides){
		LargeHashMap map = new LargeHashMap("..\\mixture_linked\\yeast_peptides.map");
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.indexPeptideByParentMass(1.0);
		map.buildTable(factory.getPeptideTable());
	}
	
	public void initialize(){
		this.factory = new TheoreticalCandidateSpectrumFactory(this.peptidesFile);
		this.queryIterator = new SortedMZXMLReader(this.queryFile);
		this.topCandidates = new TreeMap<Double, Spectrum>();
		this.pComp = SpectrumUtil.getRankBasePeakComparator(training);
		this.comp = new PRMSpectrumComparator();
		//this.comp = new SimpleProbabilisticScorer(this.pComp);
	}

	
	public void search(){
		int counter = 0;
		SimpleProbabilisticScorer comp2 = new SimpleProbabilisticScorer(this.pComp);
		comp2.setMinMatchedPeak(0);
		while(this.queryIterator.hasNext()){
			Spectrum query = this.queryIterator.next();
			query.windowFilterPeaks(8, 25);
			query.computePeakRank();
			if(query.charge > 4){
				continue;
			}
			PRMSpectrum prmq = new PRMSpectrum(query, comp2);
			ArraySpectrum aquery = new ArraySpectrum(query);
			Collection<Spectrum> candidates = factory.getCandidates(query, 3.0);
			this.topCandidates.clear();
			for(Iterator it = candidates.iterator(); it.hasNext();){
				TheoreticalSpectrum candidate = (TheoreticalSpectrum)it.next();
				double score = comp.compare(candidate, prmq);
				this.topCandidates.put(score + Math.random()*0.000000001, candidate);
			}
			printTopCandidateInfo(query, (TheoreticalSpectrum)this.topCandidates.lastEntry().getValue(), this.topCandidates.lastKey());
			counter++;
			if(counter > 1000){
				break;
			}
		}
	}
	
	private void printTopCandidateInfo(Spectrum query, TheoreticalSpectrum th, double score){
		double[] stat = th.analyzeAnnotation(query, th.peptide, 0.3);
		System.out.println("Spectrum: " + query.getPeptide() + " best-match: " + "\t" +  th.p 
			+"\t" + query.parentMass + "\t"  + query.charge + "\t" + th.parentMass + "\t" + score + "\t" 
			+ (score/th.peptide.length()) + "\t" + stat[0] + "\t" + stat[1] + "\t" + stat[2] + "\t" 
			+ stat[3] + "\t" + stat[4] + "\t" + stat[5]  + "\n");
	}
	
	public MXDBSearchSingle(){
		long start = (new GregorianCalendar()).getTimeInMillis();
		this.queryFile = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mzXML";
		this.peptidesFile = "..\\mixture_linked\\database\\Yeast_allPeptides_plusDecoy.txt";
		this.parentMassTolerance = 3.0;
		//this.fragmentMassTolerance = 0.5;
		this.initialize();
		this.search();
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void main(String[] args){
		//buildPeptideMap("..\\mixture_linked\\database\\Yeast_allPeptides_plusDecoy.txt");
		MXDBSearchSingle search = new MXDBSearchSingle();
		
	}
}
