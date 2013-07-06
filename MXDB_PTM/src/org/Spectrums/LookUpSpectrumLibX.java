package org.Spectrums;


import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A spectrum lib, indexed by peaks in the spectrum rather than
 * peptides
 * @author Jian Wang
 *
 */
public class LookUpSpectrumLibX{
	//private Map<Double, List<String>> table;
	public static int FRAGMENTMODE = 2;
	public static int PARENTMODE = 1;
	private List<String>[] table;
	private List[] pMassTable;
	List<String> pepList;
	private double massTolerance = 0.5; //default
	private static String dummy = "ABCDEFGHIJKLMNOPQQ";
	

	public LookUpSpectrumLibX(List<String> pepList){
		this(pepList, 0.5);
	}
	
	public LookUpSpectrumLibX(List<String> pepList, double massTolerance){
		this(pepList, massTolerance, FRAGMENTMODE);
//		this.pepList = pepList;
//		int maxSize = (int)((10000/massTolerance)+20);
//		this.massTolerance = massTolerance;
//		table = new ArrayList[maxSize];
//		pMassTable = new List[maxSize];
//		System.out.println("indexing " + pepList.size() + " spectra");
//		indexPeptideByPeaks(this.pepList);
//		indexPeptideByMass(this.pepList);
//		//System.out.println("after indexing we have " + table.keySet().size() + " entries");
	}
	
	public LookUpSpectrumLibX(List<String> pepList, double tolerance, int mode){
		this.pepList = pepList;
		this.massTolerance = tolerance; 
		int maxSize = (int)((10000/massTolerance)+20);
		pMassTable = new List[maxSize];
		System.out.println("indexing " + pepList.size() + " spectra");
		indexPeptideByMass(this.pepList);
		if(mode == this.FRAGMENTMODE){
			table = new ArrayList[maxSize];
			indexPeptideByPeaks(this.pepList);

		}
	}
	
	public double getMassTolerance() {
		return massTolerance;
	}

	public void setMassTolerance(double massTolerance) {
		this.massTolerance = massTolerance;
	}
	private void indexPeptideByPeaks(List<String> pepList){
		int i = 0;
		for(Iterator<String> it = pepList.iterator(); it.hasNext();){
			indexPeptide(it.next());
			i++;
			if(i % 1000 == 0){
				System.out.println("Completed " + i);
			}
			if(i%500000 == 0){
				for(int j =  0; j < this.table.length; j++){
					ArrayList current = (ArrayList)this.table[j];
					if(current != null){
						current.trimToSize();
					}
				}
			}
			if(i == 2000000){
				//break;
			}
		}
		
		for(int j =  0; j < this.table.length; j++){
			ArrayList current = (ArrayList)this.table[j];
			if(current != null){
				current.trimToSize();
			}
		}
		
		for(int j = 0; j < this.table.length; j++){
			List curr = this.table[j];
			//if(curr != null) System.out.println("entry " + j + " has size " + curr.size());
		}
	}
	
	public void indexPeptide(String peptide){
		double[][] baseMasses = TheoreticalSpectrum.computeBaseMass(peptide,	new int[]{}, new double[]{});
		for(int j = 0; j < baseMasses[0].length; j++){
				double b = baseMasses[0][j] + Mass.getIonMod("b");
				double y = baseMasses[1][j] + Mass.getIonMod("y");
				//do not index peaks with mass too small, save memory and increase efficiency
				if(b > 60)
					storeToTable(peptide, b);
				if(y > 60)
					storeToTable(peptide, y);
		}
		
	} 
	private void storeToTable(String peptide, double peakMass){
		int key=getKey(peakMass);
		int leftKey = getKey(peakMass-this.massTolerance);
		int rightKey = getKey(peakMass+this.massTolerance);
		List<String> l;
		if(this.table[key] == null){
			l = new ArrayList<String>();
		}else{
			l = this.table[key];
		}
		
		//if(!l.contains(peptide)){
			l.add(peptide);
			this.table[key] = l;
		//}
		if(leftKey != key){
			if( this.table[leftKey] == null){
				l = new ArrayList();;
			}else{
				l = this.table[leftKey];
			}
			//l.add(peptide);
			this.table[leftKey] = l;
		}
		
		if(rightKey != key && rightKey != leftKey){
			if(this.table[rightKey] == null){
				l = new ArrayList();
			}else{
				l = this.table[rightKey];
			}
			//l.add(peptide);
			this.table[rightKey] = l;
		}

	}
	
	public void indexPeptideByMass(List<String> peptides){
		for(Iterator<String> it = peptides.iterator(); it.hasNext();){
			String currentPeptide = it.next();
			indexPeptideByMass(currentPeptide);
		}
		for(int i = 0; i < this.pMassTable.length; i++){
			ArrayList curr = (ArrayList)this.pMassTable[i];
			if(curr != null){
				curr.trimToSize();
			}
		}
	}
	
	public void indexPeptideByMass(String peptide){
		//int key = getKey(PeptideMassAnalysis.computeMolecularMass(peptide));
		Peptide p = new Peptide(peptide,1);
		//System.out.println("peptide is: " + p + "\t");
		int key = getKey(p.getParentmass()-Mass.PROTON_MASS);
		List l;
		if(pMassTable[key] != null){
			l = pMassTable[key];
		}else{
			l = new ArrayList();
		}
		l.add(peptide);
		pMassTable[key] = l;
	}
	
	public int getKey(double value){
		if(value <= 0){
			return 0;
		}else{	
			//return (int)Math.round(value)+1;
			//return (int)Math.floor(value / this.massTolerance);
			return (int)Math.round((value+this.massTolerance/2)/this.massTolerance+1);
			//return (int)Math.round(value); //makes mass centered around integer
		}
	}
	
	public List<String> getSpectrumByMass(double parentMass){
		List<String> l = new ArrayList();
		if(this.pMassTable[getKey(parentMass)] != null){
			l.addAll(this.pMassTable[getKey(parentMass)]);
		}
		int key2 = getKey(parentMass-this.massTolerance);
		if(this.pMassTable[key2] != null){
			l.addAll(this.pMassTable[key2]);
		}
		int key3 = getKey(parentMass+this.massTolerance);
		if(this.pMassTable[key3] != null){
			l.addAll(this.pMassTable[key3]);
		}
		return l;
	}
	
	public List<String> getSpectrum(List<Peak> pList, int minMatch, double parentMass, double parentMassTol){
		List<String> l = new ArrayList<String>();
		Map<String, Integer> hitsCount = new HashMap<String, Integer>();
		Set<String> candidates = new HashSet();
//		double parentMass = query.parentMass*query.charge 
//			- query.charge*Mass.PROTON_MASS - Mass.WATER;
		for(int i = 0; i < pList.size(); i++){
			getSpectrumByPeak(pList.get(i), hitsCount);
		}
		for(Iterator<String> it = hitsCount.keySet().iterator(); it.hasNext();){
			String s = it.next();
			double candParentMass = PeptideMassAnalysis.computeMolecularMass(s);
			if(hitsCount.get(s).intValue() > minMatch 
					&&  Math.abs(candParentMass - parentMass) < parentMassTol){
				candidates.add(s);
			}
		}
		l.addAll(candidates);
		return l;
	}
	
	public Set<String> getSpectrumByPeak(Peak p){
		Set<String> combine = new HashSet<String>();
		return getSpectrumByPeak(p, combine);
	}
	
	
	//we should only need to lookup one key,
	//cause we tke care of error tolerance issues when indexing peaks
	//see storeToTable()

	public Set<String> getSpectrumByPeak(Peak p, Set<String> resultSet){
		int key  = getKey(p.getMass());
		//account for mass error in the real spectrum
		int key2 =getKey(p.getMass() - massTolerance);
		int key3 = getKey(p.getMass() + massTolerance);
		if(table[key] != null){
			resultSet.addAll(table[key]);
		}
//		if(key != key2 && table[key2] != null){
//			resultSet.addAll(table[key2]);
//		}
//		if(key != key3 && key2 != key3 && table[key3] != null){
//			resultSet.addAll(table[key3]);
//		}
		return resultSet;
	}
	
	public Map<String, Integer> getSpectrumByPeak(Peak p, Map<String, Integer> resultMap){
		int key  = getKey(p.getMass());
		//account for mass error in the real spectrum
		int key2 =getKey(p.getMass() - massTolerance);
		int key3 = getKey(p.getMass() + massTolerance);
		if(table[key] != null){
			incrementSpectrumCount(resultMap, table[key]);
		}
//		if(key != key2 && table[key2] != null){
//			incrementSpectrumCount(resultMap, table[key2]);
//		}
//		if(key != key3 && key2 != key3 && table[key3] != null){
//			incrementSpectrumCount(resultMap, table[key3]);
//		}
		return resultMap;
	}
	
	public Map<String, Set<Peak>> getSpectrumByPeak2(Peak p, Map<String, Set<Peak>> resultMap, Spectrum query){
		double parentMass = query.parentMass*query.charge 
		- query.charge*Mass.PROTON_MASS - Mass.WATER;
		int key  = getKey(p.getMass());
		//System.out.println("looking for peaks: " + p.getMass());
		//account for mass error in the real spectrum
		int key2 =getKey(p.getMass() - massTolerance);
		int key3 = getKey(p.getMass() + massTolerance);
		if(table[key] != null){
			//System.out.println("matched mass " + p.getMass());
			//System.out.println("getting set of size: " + table[key].size());
			incrementSpectrumCount2(resultMap, table[key], p);
		}
		
		//linked peaks
		for(int charge = getMinLinkedCharge(query.charge); charge <= getMaxLinkedCharge(query.charge); charge++){
			//double complement = (parentMass - (p.getMass()*charge-Mass.PROTON_MASS*(charge-1))) + Mass.getIonMod("b") + Mass.getIonMod("y");
			double complement = (parentMass - (p.getMass()*charge-Mass.PROTON_MASS*(charge-1))) + Mass.getIonMod("b") + Mass.getIonMod("y");			
			//System.out.println("complement is: " + complement + " peak mass is: " + p.getMass());
			int keyComplement = getKey(complement);
			//System.out.println("looking for complement peaks: " + complement);
			//we need to expand mass error tolerance for x-linked peptide 
			//since complement fragment can have error as large as  tolerance*peakcharge
			for(int i = -1; i <= 1; i++){
				int keyComplement2 = getKey(complement+i*this.massTolerance);
				if(table[keyComplement2] != null){
					//incrementSpectrumCount2(resultMap, table[keyComplement2], p);
				}
			}
		}
		return resultMap;
	}
	
	private int getMinLinkedCharge(int linkedCharge){
		if(linkedCharge <= 3){
			return 2;
		}else{
			return 3;
		}
	}
	
	private int getMaxLinkedCharge(int linkedCharge){
		return linkedCharge-1;
	}
	
	public Map<Spectrum, Integer> getSpectrumByLinkedPeak(Peak p, Map<Spectrum, Integer> resultMap, double parentMass){
		double mass = p.getMass();
		double mass2 = p.getMass()*2;
		//double 
		return resultMap;
	}
	
	private void incrementSpectrumCount(Map<String, Integer> m, List<String> s){
		for(Iterator<String> it = s.iterator(); it.hasNext();){
			String curr = it.next();
			if(m.containsKey(curr)){
				Integer count = m.get(curr);
				m.put(curr, count+1);
			}else{
				m.put(curr, new Integer(1));
			}
		}
	}
	
	private void incrementSpectrumCount2(Map<String, Set<Peak>> m, List<String> s, Peak p){
		for(Iterator<String> it = s.iterator(); it.hasNext();){
			String curr = it.next();
			if(m.containsKey(curr)){
				Set<Peak> matched = m.get(curr);
				matched.add(p);
			}else{
				Set<Peak> matched = new HashSet();
				matched.add(p);
				m.put(curr, matched);
			}
		}
	}
	
	public List<String> getSpectrumByPeaks(List<Peak> peakList){
		Set<String> combined = new HashSet<String>();
		List<String> l = new ArrayList<String>();
		for(int i = 0; i < peakList.size(); i++){
			getSpectrumByPeak(peakList.get(i), combined);
		}
		l.addAll(combined);
		return l;
	}
	
	//for regular peptides
	public List<String> getSpectrumByPeaks(List<Peak> peakList, int minMatch){
		List<String> l = new ArrayList<String>();
		Map<String, Integer> hitsCount = new HashMap<String, Integer>();
		Set<String> candidates = new HashSet();
//		double parentMass = query.parentMass*query.charge 
//			- query.charge*Mass.PROTON_MASS - Mass.WATER;
		for(int i = 0; i < peakList.size(); i++){
			getSpectrumByPeak(peakList.get(i), hitsCount);
		}
		for(Iterator<String> it = hitsCount.keySet().iterator(); it.hasNext();){
			String s = it.next();
			if(hitsCount.get(s).intValue() > minMatch){
				candidates.add(s);
			}
		}
		l.addAll(candidates);
		return l;
	}
	
	//for linked peptides
	public List<String> getSpectrumByPeaks(List<Peak> peakList, int minMatch, Spectrum query){
		List<String> l = new ArrayList<String>();
		//Map<String, Integer> hitsCount = new HashMap<String, Integer>();
		Map<String, Set<Peak>> hitsCount = new HashMap<String, Set<Peak>>();
		Set<String> candidates = new HashSet();
		double parentMass = query.parentMass*query.charge 
			- query.charge*Mass.PROTON_MASS - Mass.WATER;
		for(int i = 0; i < peakList.size(); i++){
			getSpectrumByPeak2(peakList.get(i), hitsCount, query);
		}
		for(Iterator<String> it = hitsCount.keySet().iterator(); it.hasNext();){
			String s = it.next();
			if(hitsCount.get(s).size() > minMatch){
				//getLinkedPair(s, minMatch, query, hitsCount, candidates);
				candidates.add(s);
			}
		}
		l.addAll(candidates);
		return l;
	}
	
	
	private void getLinkedPair(String peptide, int minMatch, Spectrum query, Map<String, Set<Peak>> hitsCount, Set<String> cands){
		double parentMass = query.parentMass * query.charge - query.charge*Mass.PROTON_MASS; //deionized masses
		double thisMass = PeptideMassAnalysis.computeMolecularMass(peptide);
		double partnerMass = parentMass - thisMass - Mass.DSPLINKER_MASS;
		//System.out.println("query : " + query.peptide + " one cand: " + peptide + " parentmass: " + parentMass + " this mass " + thisMass + "partnermass: " + partnerMass);
		double tolerance = 0.5;
		int key = getKey(partnerMass);
		int pairCandCount = 0;
		if(this.pMassTable[key] != null){
			List<String> paircand = this.pMassTable[key];
			for(Iterator<String> it = paircand.iterator(); it.hasNext();){
				String currPair = it.next();
				if(hitsCount.containsKey(currPair) && 
						hitsCount.get(currPair).size() + hitsCount.get(peptide).size() > 0){
					double candmass = PeptideMassAnalysis.computeMolecularMass(currPair);
					//System.out.println("query : " + query.peptide + " cand: " + peptide + " and " + currPair + 
					//		" has presumed mass: " + (candmass + thisMass + Mass.DSPDANGLE_MASS) + " parentmass: " + parentMass);
					if(Math.abs(candmass + thisMass + Mass.DSPLINKER_MASS - parentMass) < tolerance){
						cands.add(currPair);
						pairCandCount++;
					}
				}
			}
			if(pairCandCount > 1){
				cands.add(peptide);
			}
		}
	}
	
	public static double getLinkedOffSet(String peptide, Spectrum linkedquery){
		double parentMass = linkedquery.parentMass * linkedquery.charge - linkedquery.charge*Mass.PROTON_MASS; //deionized masses
		double thisMass = PeptideMassAnalysis.computeMolecularMass(peptide);
		double partnerMass = parentMass - thisMass;
		return partnerMass;
	}
	
	public static double getLinkedOffSet(Peptide peptide, Spectrum linkedquery){
		double parentMass = linkedquery.parentMass * linkedquery.charge - linkedquery.charge*Mass.PROTON_MASS; //deionized masses
		double thisMass = (peptide.getParentmass()*peptide.getCharge() - peptide.getCharge()*Mass.PROTON_MASS);
		double partnerMass = parentMass - thisMass;
		return partnerMass;
	}
	
	//need to check whether crosslinker has correct mass implementation
	public static double getLinkedPartnerParentmass(String peptide, Spectrum linkedquery, CrossLinker l){
		return getLinkedOffSet(peptide, linkedquery) - l.getLinkerMassOffSet();
		
	}
	
	public static double getLinkedPartnerParentmass(Peptide peptide, Spectrum linkedquery, CrossLinker l){
		for(int i = 0; i < peptide.getPos().length; i++){
			if(peptide.getPos()[i] == peptide.getLinkedPos()){
				return peptide.getPtmmasses()[i]-l.getLinkerMassOffSet();
			}
		}
		return -1;
	}
	
	public static boolean checkPassFilter(String peptide, List<String> filtered){
		for(int i = 0; i < filtered.size(); i++){
			if(peptide.equals(filtered.get(i))){
				return true;
			}
		}
		return false;
	}
	
	public static int checkPassFilter(String peptide1, String peptide2, List<String> filtered){
		int matchCount = 0;
		System.out.println("peptides are: " + peptide1 + " " + peptide2);
		for(int i = 0; i < filtered.size(); i++){
			String curr = filtered.get(i);
			if(curr.equals(peptide1)){
				matchCount+=1;
			}
			if(curr.equals(peptide2)){
				matchCount+=2;
			}
		}
		return matchCount;
	}
	
	public static void testLookUpSpectrumX(){
		String spectrumFile = "..\\mixture_linked\\linked_peptide_spectra1.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
//		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
//		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		System.out.println("done loading query spectra lib");
		String file = "..\\mixture_linked\\database\\Yeast_allPeptides_plusLinked_peptides_plusDecoy.txt";
		//String file = "..\\mixture_linked\\tempLinkedPeptides.txt";
		lib1.removeModSpectra();
		List<Spectrum> specList = lib1.getSpectrumList();
		List<String> pepList = Utils.FileIOUtils.createListFromFile(file);
		LookUpSpectrumLibX lookup = new LookUpSpectrumLibX(pepList);
//		factory = null;
		System.out.println("Done indexing peptides");
		long start = (new GregorianCalendar()).getTimeInMillis();
		//SpectrumLib lib = factory.createCandidateSpectrumLibX(lib1.getRandomSpectrum(), 10000, false);
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			s.windowFilterPeaks(10, 25);
			//s.computePeakRank();
			List<Peak> queryPeaks = s.getTopPeaks(25);
			//List<Peak> queryPeaks = s.topWindowPeaks(6, 200);
			System.out.println("Query peaks has: " + queryPeaks.size());
			List<String> candidates = lookup.getSpectrumByPeaks(queryPeaks, 3, s);
			//List<String> candidates = lookup.getSpectrumByPeaks(queryPeaks, 2);
//			String[] peptides = s.peptide.split(" & ");
//			int passedFilter = checkPassFilter(peptides[0], peptides[1].substring(0, peptides[1].length()-2), candidates);
			System.out.println("Query: "  + s.spectrumName + " After filter one we have: " + candidates.size() + " candidates ");
//			System.out.println("After filter correct peptide is retained?: " + passedFilter);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testIndexingMemory(){
		List[] table = new List[7000];
		for(int i = 0; i < table.length; i++){
			table[i] = new ArrayList(10000);
		}
		System.out.println("start adding");
		for(int i = 0; i < table.length; i++){
			List<String> pepList = table[i];
			for(int j = 0; j < 10000; j++){
				pepList.add(null);
			}
		}
		System.out.println("Done adding");
	}
	
	public static void testIndexingMemory2(){
		List<Double> peptides = new ArrayList<Double>(70000000);
		for(int i = 0; i < 70000000; i++){
			peptides.add(Math.random());
		}
	}
	
	public static void main(String[] args){
		testLookUpSpectrumX();
		//testIndexingMemory();
		//testIndexingMemory2();
	}
}
