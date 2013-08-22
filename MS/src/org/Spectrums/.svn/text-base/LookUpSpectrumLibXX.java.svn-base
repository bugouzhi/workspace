package org.Spectrums;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * Algorithm base on paper by Nathan Edward and Ross  Lippert LNCS 2452, pp. 68-81, 2002
 * Use to generate peptide candidates that match precursor mass and at least K b/y ions
 * to the top N most intense peaks in the database
 * @author bugouzhi
 *
 */
public class LookUpSpectrumLibXX {
	private String proteins;
	private String proteinFileName;
	private int[] siteIndex;
	private int toleranceMode = 0; //0: Da; 1: ppm
	private double parentMassTolerance = 3000.0;
	private double fragmentMassTolerance = 1.0;
	private int minMatchedPeak = 3;
	private int minContinuousMatch = 1;
	private int minCharge = 2;
	private int maxCharge = 4;
	private int minLinkedCharge = 2;
	private int maxLinkedCharge = 6;
	private Set<Peak> exclusion = new HashSet<Peak>();
	private boolean DEBUG = false;
    public double[] modList= new double[]{0.0}; //always put +0.0 mod to stand for non-modified peptides
    public int maxMods=1;
    public int numC13=1;
	
	public LookUpSpectrumLibXX(){
    	
    }
    
	public int getToleranceMode() {
		return toleranceMode;
	}

	public void setToleranceMode(int toleranceMode) {
		this.toleranceMode = toleranceMode;
	}
	
    public int getMinMatchedPeak() {
		return minMatchedPeak;
	}

	public void setMinMatchedPeak(int minMatchedPeak) {
		this.minMatchedPeak = minMatchedPeak;
	}

	public int getMinContinuousMatch() {
		return minContinuousMatch;
	}

	public void setMinContinuousMatch(int minContinuousMatch) {
		this.minContinuousMatch = minContinuousMatch;
	}


	public int getMinCharge() {
		return minCharge;
	}

	public void setMinCharge(int minCharge) {
		this.minCharge = minCharge;
	}

	public int getMaxCharge() {
		return maxCharge;
	}

	public void setMaxCharge(int maxCharge) {
		this.maxCharge = maxCharge;
	}

	public int getMinLinkedCharge() {
		return minLinkedCharge;
	}

	public void setMinLinkedCharge(int minLinkedCharge) {
		this.minLinkedCharge = minLinkedCharge;
	}

	public int getMaxLinkedCharge() {
		return maxLinkedCharge;
	}

	public void setMaxLinkedCharge(int maxLinkedCharge) {
		this.maxLinkedCharge = maxLinkedCharge;
	}

	public static void main(String[] args){
    	//testLoadProteins();
    	testFiltering();
    }
    //concat peptides into one long string, record index
    //as if cut by enzyme
    public void loadPeptidesFromFile(String peptideFile){
    	List<Integer> index = new ArrayList<Integer>();
    	index.add(0);
    	StringBuffer buff = new StringBuffer();
    	try{
    		BufferedReader reader = new BufferedReader(new FileReader(peptideFile));
    		String line = reader.readLine();
    		while(line != null){
    			buff.append(line);
    			buff.append("*");
    			index.add(buff.length()-1);
    			line = reader.readLine();
    		}
    		buff.append("*");//append dummy residue at the end
    		index.add(buff.length()-1);
    		this.siteIndex = new int[index.size()];
    		for(int i = 0; i < index.size(); i++){
    			siteIndex[i] = index.get(i);
    		}
    		System.out.println("Loaded peptides: " + index.size());
    		this.proteinFileName = peptideFile;
    		this.proteins = buff.toString();
    	}catch(IOException ioe){
    		System.out.println("Error reading file");
    		System.out.println(ioe.getMessage());
    		ioe.getStackTrace();
    	}
    }
    
    public void loadPeptidesFromFileLite(String peptideFile){
    	List<Integer> index = new ArrayList<Integer>();
    	index.add(0);
    	StringBuffer buff = new StringBuffer();
    	try{
    		BufferedReader reader = new BufferedReader(new FileReader(peptideFile));
    		String line = reader.readLine();
    		while(line != null){
    			buff.append(line);
    			line = reader.readLine();
    		}
    		for(int i = 0; i < buff.length(); i++){
    			if(buff.charAt(i) == '*'){
    				index.add(i);
    			}
    		}
    		this.siteIndex = new int[index.size()];
    		for(int i = 0; i < index.size(); i++){
    			siteIndex[i] = index.get(i);
    		}
    		System.out.println("Loaded peptides: " + index.size());
    		this.proteinFileName = peptideFile;
    		this.proteins = buff.toString();
    	}catch(IOException ioe){
    		System.out.println("Error reading file");
    		System.out.println(ioe.getMessage());
    		ioe.getStackTrace();
    	}
    }
    
    public List<PeptideLite> getCandidatePeptide(double parentMass, int charge, List<Peak> pList){
    	double currentMass = 0;
    	int currentSite = 1;
    	int currentIndex = 0;
    	int totalLength = proteins.length();
    	//List<Peak> complements = this.getLinkedComplementPeaks(pList, parentMass, charge);
    	//System.out.println("number of filtering peaks " + pList.size());
    	Map<Integer, Peak> pLookUpTable = createMatchedPeakTable(pList, 1);
    	createComplementPeakTable(pList, parentMass, charge, 2, pLookUpTable);
    	this.exclusion.clear();
    	int matchedPeakCount = 0;
    	int continuousMatch=0;
    	List<PeptideLite> candidates = new ArrayList<PeptideLite>();
    	if(parentMass <= 71){
    		return candidates;
    	}
    	//System.out.println("parentmass: " + parentMass);
    	while(currentIndex < totalLength && currentSite < siteIndex.length){
    		while((currentMass/charge) <= (parentMass + this.parentMassTolerance) 
    				&& currentIndex < this.siteIndex[currentSite]){
    			currentMass += Mass.getAAMass(proteins.charAt(currentIndex));
    			if(DEBUG){
    				System.out.println("current mass: " + currentMass  + "\t"
    					+ this.proteins.substring(this.siteIndex[currentSite-1]+1, siteIndex[currentSite]));
    			}
//    			if(checkPeakLookUpTable(pLookUpTable, currentMass+Mass.getIonMod("b"))){ 
//    					//|| checkPeakLookUpTable(pLookUpTable, currentMass+Mass.getIonMod("b-NH3"))){
//    				if(DEBUG){
//    					System.out.println("matched peak");
//    				}
//    				matchedPeakCount++;
//    				continuousMatch++;
//    			}else{
//    				if(continuousMatch <= this.minContinuousMatch){
//    					continuousMatch = 0;
//    				}
//    			}
    			//System.out.println("currentIndex " + currentIndex);
    			//System.out.println("currentMass " + currentMass);
    			currentIndex++;
    		} 
    		if(continuousMatch <= this.minContinuousMatch){
    			continuousMatch = 0;
    		}
    		if(currentIndex == siteIndex[currentSite]){
    			//System.out.println("checking parentmasss: " + parentMass + "\t" + currentMass);
//    			if(this.proteins.substring(siteIndex[currentSite-1]+1, 
//    								siteIndex[currentSite]).startsWith("KAADAE")){
//    				System.out.println("checking: " + this.proteins.substring(siteIndex[currentSite-1]+1, 
//    								siteIndex[currentSite]));
//    			}
    			if(checkParentMass(parentMass, currentMass, charge)){
    				//System.out.println("matching parentmass");
    				if(matchedPeakCount > this.minMatchedPeak && 
    							continuousMatch > this.minContinuousMatch){
    					candidates.add(new PeptideLite(siteIndex[currentSite-1]+1, 
    								siteIndex[currentSite],
    								this.proteins, 
    								2));
    					//System.out.println("adding: " + this.proteins.substring(siteIndex[currentSite-1]+1,siteIndex[currentSite]));
    				}else{  //we travel backward to match y-ions
    					int backIndex = currentIndex-1;
    					double backWardMass = 0;
    					while(backIndex >= siteIndex[currentSite-1]){
    						backWardMass += Mass.getAAMass(proteins.charAt(backIndex));
    						if(DEBUG){
    							System.out.println("current mass: " + backWardMass  + "\t"
    									+ this.proteins.substring(backIndex+1, siteIndex[currentSite]));
    						}
    						if(checkPeakLookUpTable(pLookUpTable, backWardMass+Mass.getIonMod("y"))){
    								//|| checkPeakLookUpTable(pLookUpTable, backWardMass+Mass.getIonMod("y-NH3"))){
    							if(DEBUG){
    								System.out.println("matched peaks");
    							}
    		    				matchedPeakCount++;
    		    				continuousMatch++;
    		    			}else{
    		    				if(continuousMatch<=this.minContinuousMatch){
    		    					continuousMatch = 0;
    		    				}
    		    			}
    						if(matchedPeakCount > this.minMatchedPeak 
    								&& continuousMatch > this.minContinuousMatch){
    							candidates.add(new PeptideLite(siteIndex[currentSite-1]+1, 
        								siteIndex[currentSite],
        								this.proteins,
        								2));
    							break;
    						}
    						//System.out.println("currentback: " + backIndex);
    						backIndex--;
    					}
    				}
    			}
    		}
    		currentIndex = siteIndex[currentSite]+1;
			currentSite++;
			//System.out.println("currentSite: " + currentSite);
			currentMass = 0;
			matchedPeakCount=0;
			continuousMatch =0;
			this.exclusion.clear();
    	}
    	return candidates;
    }
    
    public double getParentMassTolerance() {
		return parentMassTolerance;
	}

	public void setParentMassTolerance(double parentMassTolerance) {
		this.parentMassTolerance = parentMassTolerance;
	}

	public double getFragmentMassTolerance() {
		return fragmentMassTolerance;
	}

	public void setFragmentMassTolerance(double fragmentMassTolerance) {
		this.fragmentMassTolerance = fragmentMassTolerance;
	}
    
	private boolean checkPeakLookUpTable(Map<Integer, Peak> lookup, double mass){
    	if(lookup.containsKey(getKey(mass))){
    		Peak value = lookup.get(getKey(mass));
    		if(!this.exclusion.contains(value)){
    			this.exclusion.add(value);
    			return true;
    		}
    	}
    	return false;
    }
	
    private boolean checkParentMass(double parentMass, double currentMass, int charge){
    	for(int c = this.minCharge; c <= this.maxCharge; c++){
    		double candidateMass = (currentMass + Mass.WATER + Mass.PROTON_MASS*charge)/charge;
    		for(int i = 0; i < this.modList.length; i++){
    			for(int j =0; j <= this.numC13; j++){
    				candidateMass = (currentMass + Mass.WATER + modList[i] + (Mass.C13-Mass.C12)*j + Mass.PROTON_MASS*c)/c;
    				//System.out.println("error: " + (parentMass - candidateMass));
    				if(this.toleranceMode == 0){
    					if(Math.abs(parentMass - candidateMass) < this.parentMassTolerance*charge){
    						return true;
    					}
    				}else{
    					//System.out.println(parentMass + "\t" + candidateMass + "\t" + "error: " + (parentMass - candidateMass)*1000000/parentMass);
    					if(Math.abs(parentMass - candidateMass)*1000000/parentMass < this.parentMassTolerance){
    						return true;
    					}
    				}
    			}
    		}
    	}
    	return false;
    }
    
    private Map<Integer, Peak> createMatchedPeakTable(List<Peak> pList, int width){
    	Map<Integer, Peak> pLookUpTable = new HashMap<Integer, Peak>();
    	return createMatchedPeakTable(pList, width, pLookUpTable);
    }
    
    private Map<Integer, Peak> createMatchedPeakTable(List<Peak> pList, int width, Map<Integer,Peak> pLookUpTable){
    	for(Iterator<Peak> it = pList.iterator(); it.hasNext();){
    		Peak current = it.next();
    		pLookUpTable.put(getKey(current.getMass()), current);
    		for(int i = 1; i < width; i++){
    			pLookUpTable.put(getKey(current.getMass()-i*this.fragmentMassTolerance), current);
    			pLookUpTable.put(getKey(current.getMass()+i*this.fragmentMassTolerance), current);
    		}
    	}
    	return pLookUpTable;
    }
    
    private Map<Integer, Peak> createComplementPeakTable(List<Peak> pList, double parentMass, 
    			int linkedCharge, int width, Map<Integer,Peak> pLookUpTable){
    	for(Iterator<Peak> it = pList.iterator(); it.hasNext();){
    		Peak current = it.next();
    		getLinkedComplementPeaks(current, parentMass, linkedCharge, width, pLookUpTable);
    	}
    	return pLookUpTable;
    }
    
    //assume around integer values, may not suitable for high mass accuracy data
	public int getKey(double value){
		//value = value*0.9995; //center masses values around integer, add in later
		if(value <= 0){
			return 0;
		}else{	
			//return (int)Math.round(value)+1;
			//return (int)Math.floor(value / this.massTolerance);
			return (int)Math.round((value+this.fragmentMassTolerance/2)/this.fragmentMassTolerance+1);
			//return (int)Math.round(value); //makes mass centered around integer
		}
	}
	
	public List<Peak> getLinkedComplementPeaks(List<Peak> pList, double parentMass, int linkedCharge){
		parentMass =  parentMass*linkedCharge 
			- linkedCharge*Mass.PROTON_MASS - Mass.WATER;
		if(DEBUG){
			System.out.println("parentmass: " + parentMass);
		}
		int length = pList.size();
		List<Peak> complements = new ArrayList<Peak>();
		for(int i = 0; i < length; i++){
			Peak p = pList.get(i);
			//linked peaks
			for(int charge = LinkedPeptide.getMinLinkedCharge(linkedCharge); 
				charge <= LinkedPeptide.getMaxLinkedCharge(linkedCharge); charge++){
				//double complement = (parentMass - (p.getMass()*charge-Mass.PROTON_MASS*(charge-1))) + Mass.getIonMod("b") + Mass.getIonMod("y");
				double complement = (parentMass - (p.getMass()*charge-Mass.PROTON_MASS*(charge-1))) + Mass.getIonMod("b") + Mass.getIonMod("y");
				if(true){
					System.out.println("storing complement peaks: " + complement + " original peaks: " + p);
				}
				complements.add(new Peak(complement, p.getIntensity()));
			}
		}
		return complements;
	}
	
	public void getLinkedComplementPeaks(Peak p, double parentMass, int linkedCharge, int width, Map<Integer,Peak> pLookUpTable){
		parentMass =  parentMass*linkedCharge 
			- linkedCharge*Mass.PROTON_MASS - Mass.WATER; 
		//linked peaks
		for(int charge = LinkedPeptide.getMinLinkedCharge(linkedCharge); 
			charge <= LinkedPeptide.getMaxLinkedCharge(linkedCharge); charge++){
			//double complement = (parentMass - (p.getMass()*charge-Mass.PROTON_MASS*(charge-1))) + Mass.getIonMod("b") + Mass.getIonMod("y");
			double complement = (parentMass - (p.getMass()*charge-Mass.PROTON_MASS*(charge-1))) + Mass.getIonMod("b") + Mass.getIonMod("y");
			double complement2 = (parentMass - (p.getMass()*charge-Mass.PROTON_MASS*(charge-1))) + Mass.getIonMod("b") + Mass.getIonMod("y") - Mass.NH3;
			pLookUpTable.put(getKey(complement), p);
			//pLookUpTable.put(getKey(complement2), p);
			if(DEBUG){
				System.out.println("storing complement peaks: " + complement + " original peaks: " + p);
				System.out.println("storing complement peaks2: " + complement2 + " original peaks: " + p);
			}
			for(int i = 1; i < width; i++){
    			pLookUpTable.put(getKey(complement-i*this.fragmentMassTolerance), p);
    			pLookUpTable.put(getKey(complement+i*this.fragmentMassTolerance), p);
    			pLookUpTable.put(getKey(complement2-i*this.fragmentMassTolerance), p);
    			pLookUpTable.put(getKey(complement2+i*this.fragmentMassTolerance), p);

    		}
			
		}
	}
	
	public static int checkPassFilter(String peptide1, String peptide2, List<PeptideLite> filtered){
		int matchCount = 0;
		boolean found1 = false, found2 = false;
		System.out.println("peptides are: " + peptide1 + " " + peptide2);
		for(int i = 0; i < filtered.size(); i++){
			String curr = filtered.get(i).getPep();
			//System.out.println("checking peptide " + filtered.get(i).getPep());
			if(curr.equals(peptide1)){
				found1 = true;
			}
			if(curr.equals(peptide2)){
				found2 = true;
			}
		}
		if(found1 && found2){
			return 3;
		}else if(found1){
			return 1;
		}else if(found2){
			return 2;
		}else{
			return 0;
		}
	}
    public static void testLoadProteins(){
    	String peptideFile = "..\\mixture_linked\\database\\Human_allpeptides_plusLinkedPeptdies_plusDecoy.txt";
    	String spectrumFile = "..\\mixture_linked\\linked_peptide_spectra2.mgf";
    	SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
    	LookUpSpectrumLibXX lookup = new LookUpSpectrumLibXX();
    	lookup.loadPeptidesFromFile(peptideFile);
    	List<Spectrum> specList = lib1.getAllSpectrums();
    	System.out.println("start searching");
    	long start = (new GregorianCalendar()).getTimeInMillis();
    	for(int i = 0; i < specList.size(); i++){
    		Spectrum s = specList.get(i);
    		s.windowFilterPeaks(10, 25);
    		List<Peak> pList = s.getTopPeaks(20);
    		System.out.println("Searching " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
    		List<PeptideLite> candidates = lookup.getCandidatePeptide(s.parentMass, s.charge, pList);
    		String[] peptides = s.peptide.split(" & ");
			int passedFilter = checkPassFilter(peptides[0], peptides[1], candidates);
    		System.out.println(s.spectrumName + " has candidates: " + candidates.size());
			System.out.println("After filter correct peptide is retained?: " + passedFilter);
 
    	}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
    }
    
    
    public static void testFiltering(){
    	String peptideFile = "..\\mixture_linked\\database\\Human_allpeptides_plusLinkedPeptdies_plusDecoy.txt";
    	String spectrumFile = "..\\mixture_linked\\linked_peptide_spectra1.mgf";
    	String libfile2 = "..\\MSPLib\\Lib\\ecoli.msp";
    	SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
    	LookUpSpectrumLibXX lookup = new LookUpSpectrumLibXX();
    	lookup.loadPeptidesFromFile(peptideFile);
    	List<Spectrum> specList = lib1.getAllSpectrums();
    	System.out.println("start searching");
    	long start = (new GregorianCalendar()).getTimeInMillis();
    	SpectrumComparator scorer = SpectrumUtil.getLPeakRankBaseScorer(libfile2);
    	for(int i = 0; i < specList.size(); i++){
    		Spectrum s = specList.get(i);
    		s.windowFilterPeaks(10, 25);
    		s.computePeakRank();
    		List<Peak> pList = s.getTopPeaks(20, 1.0);
    		System.out.println("Searching " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
    		List<PeptideLite> candidates = lookup.getCandidatePeptide(s.parentMass, s.charge, pList);
    		String[] peptides = s.peptide.split(" & ");
			int passedFilter = checkPassFilter(peptides[0], peptides[1], candidates);
    		System.out.println(s.spectrumName + " has candidates: " + candidates.size());
			System.out.println("After filter correct peptide is retained?: " + passedFilter);
			List<Spectrum> candidateSpectrum = LinkedPeakScoreLearner.generateSpectra2(candidates, s);
			for(int j = 0; j < candidateSpectrum.size(); j++){
				((TheoreticalSpectrum)candidateSpectrum.get(j)).getPeak();
			}
			SpectrumLibSearcher searcher = new SpectrumLibSearcher(candidateSpectrum, scorer);
			int[] ranks = searcher.linkedRanks(s);
			System.out.println("target peptides ranks " + ranks[0] + "\t" + ranks[1]);
    	}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
    }

}
