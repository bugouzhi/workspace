package org.Spectrums;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Generate theoretical spectrum candidates;
 * back by a Peptide candidates table stored in file
 * @author Jian Wang
 *
 */
public class TheoreticalCandidateSpectrumFactory {
	LargeHashMap peptideMap;
	TreeSet<Spectrum> candidates;
	String peptidesFile;
	private double currentMax;
	private boolean matchCharge = false;
	public boolean isMatchCharge() {
		return matchCharge;
	}

	public void setMatchCharge(boolean matchCharge) {
		this.matchCharge = matchCharge;
	}

	public TheoreticalCandidateSpectrumFactory(CandidateSpectrumLibFactory factory){
		this.peptidesFile = factory.peptideFile +"_processed.map";
		System.out.println("peptide file is: " + this.peptidesFile);
		System.out.println("factory has " + factory.peptideTable.values().size() + " peptides");
		LargeHashMap map = new LargeHashMap(this.peptidesFile);
		map.buildTable(factory.getPeptideTable());
		init();
	}
	
	public TheoreticalCandidateSpectrumFactory(String peptidesFile){
		this.peptidesFile = peptidesFile;
		if(!peptidesFile.endsWith(".map")){
			this.peptidesFile = peptidesFile+"_processed.map";
			LargeHashMap map = new LargeHashMap(peptidesFile+"_processed.map");
			CandidateSpectrumLibFactory factory = 	
				CandidateSpectrumLibFactory.createFactoryFromPeptide(peptidesFile);
			factory.indexPeptideByParentMass(1.0);
			map.buildTable(factory.getPeptideTable());
		}
		init();
	}
	
	private void init(){
		this.peptideMap = new LargeHashMap(this.peptidesFile);
		this.peptideMap.loadLibraryFromFile(this.peptidesFile);
		this.candidates = new TreeSet(SpectrumMassComparator.comparator);
		this.currentMax = 0.0;
	}
	
	
	public Collection<Spectrum> getCandidates(Spectrum s, double parentMassTolerance){
		long key = getKey(s);
		long width = (long)Math.ceil(parentMassTolerance);
		System.out.println("key is : " + key + "\twidth:\t" + width);
		TreeSet<Peptide> candPeptides = new TreeSet<Peptide>(PeptideMassComparator.comparator);
		for(long i = key-width; i <= key+width; i++){
			List<Peptide> value = (List<Peptide>)this.peptideMap.get(i);
			if(value != null){
				candPeptides.addAll(value);
			}
		}
		//remove candidates
		Spectrum current=null;
		System.out.println("number of candidate before: " + this.candidates.size());
		for(Iterator iter = this.candidates.iterator(); iter.hasNext();){
			current = (Spectrum)iter.next();
			if(Math.abs(s.parentMass - current.parentMass) > parentMassTolerance 
					|| (matchCharge && s.charge != current.charge)){
				iter.remove();
			}
		}
		System.out.println("number of candidate after: " + this.candidates.size());
		System.out.println("number of candidate peptides: " + candPeptides.size());
		Peptide p = new Peptide();
		if(current == null){
			p.setParentmass(0.0000001);
		}else{
			p.setParentmass(current.parentMass);
		}
		System.out.println("current top key: " + p.getParentmass());
		System.out.println("adding trailing candidates: " + candPeptides.tailSet(p).size());
		//this.candidates.clear();
		int added=0;
		for(Iterator iter = candPeptides.tailSet(p).iterator(); iter.hasNext();){
			Peptide cand = (Peptide)iter.next();
			if(Math.abs(cand.getParentmass() - s.parentMass) <= parentMassTolerance 
					&& (!matchCharge || s.charge == current.charge)){
				this.candidates.add(new ArraySpectrum(new TheoreticalSpectrum(cand)));
				//this.candidates.add(new TheoreticalSpectrum(cand));
				added++;
			}	
		}
		System.out.println("added candidates: " + added);
		System.out.println("candidates has size: " + this.candidates.size());
		//List<TheoreticalSpectrum> l = new ArrayList();
		//l.addAll(this.candidates);
		return this.candidates;
	}
	
	public Collection<Spectrum> getSUMOCandidates(Spectrum s, double parentMassTolerance){
		List<Spectrum> cand = new ArrayList();
		cand.addAll(getSUMOCandidates(s.parentMass, s.charge, parentMassTolerance));
		//cand.addAll(getSUMOCandidates(s.parentMass+1.0/s.charge, s.charge, parentMassTolerance));
		cand.addAll(getSUMOCandidates(s.parentMass-1.0/s.charge, s.charge, parentMassTolerance));
		//return getSUMOCandidates(s.parentMass, s.charge, parentMassTolerance);
		return cand;
	}
	public Collection<Spectrum> getSUMOCandidates(double parentMass, int charge, double parentMassTolerance){
		long key = getKey(parentMass);
		long width = (long)Math.ceil(parentMassTolerance);
		System.out.println("key is : " + key + "\twidth:\t" + width);
		TreeSet<Peptide> candPeptides = new TreeSet<Peptide>(PeptideMassComparator.comparator);
		for(long i = key-width; i <= key+width; i++){
			List<Peptide> value = (List<Peptide>)this.peptideMap.get(i);
			if(value != null){
				candPeptides.addAll(value);
			}
		}
		//remove candidates
		Spectrum current=null;
		System.out.println("number of candidate before: " + this.candidates.size());
		for(Iterator iter = this.candidates.iterator(); iter.hasNext();){
			current = (Spectrum)iter.next();
			if(Math.abs(parentMass - current.parentMass) > parentMassTolerance){
				iter.remove();
			}
		}
		System.out.println("number of candidate after: " + this.candidates.size());
		System.out.println("number of candidate peptides: " + candPeptides.size());
		//System.out.println("current is: " + current);
		Peptide p = new Peptide();
		if(current == null){
			p.setParentmass(0.0000001);
		}else{
			p.setParentmass(current.parentMass);
		}
		System.out.println("current top key: " + p.getParentmass() + "\t" + p.getPeptide());
		System.out.println("adding trailing candidates: " + candPeptides.tailSet(p).size());
		//this.candidates.clear();
		int added=0;
		for(Iterator iter = candPeptides.tailSet(p).iterator(); iter.hasNext();){
			Peptide cand = (Peptide)iter.next();
			if(Math.abs(cand.getParentmass() - parentMass) <= parentMassTolerance){
				//TheoreticalSpectrum t = new TheoreticalSpectrum((LinkedPeptide)cand, (short)cand.getCharge(), false);
				LazyEvaluateLinkedSpectrum t = new LazyEvaluateLinkedSpectrum(((LinkedPeptide)cand).peptides[0], 
						((LinkedPeptide)cand).peptides[1], cand.getCharge(), Mass.DSSLINKER_MASS);
				//System.out.println(cand);
				//System.out.println(t);
				this.candidates.add(t);
				//this.candidates.add(new TheoreticalSpectrum(cand));
				added++;
			}	
		}
		System.out.println("added candidates: " + added);
		System.out.println("candidates has size: " + this.candidates.size());
		Collection<Spectrum> l = new ArrayList();
		for(Iterator<Spectrum> it = this.candidates.iterator(); it.hasNext();){
			Spectrum spect = it.next();
			if(!this.matchCharge || spect.charge == charge)
				l.add(spect);
		}	
		//l.addAll(this.candidates);
		return l;
	}
	
	public List<Peptide> getPeptideByMass(double mass){
		int key = getKey(mass);
		List<Peptide> cand = (List<Peptide>)this.peptideMap.get(key);
		return cand;
	}
	
	private long getKey(Spectrum s){
		return (long)Math.round(s.parentMass/1.0);
	}
	
	private int getKey(double parentMass){
		return (int)Math.round(parentMass/1.0);
	}
	
	public static void buildCandidateTable(CandidateSpectrumLibFactory factory, String file){
		LargeHashMap map = new LargeHashMap(file);
		factory.indexPeptideByParentMass(1.0);
		map.buildTable(factory.getPeptideTable());
	}
	
	
	
}
