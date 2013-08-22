package org.Spectrums;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
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
 * Generate candidates like CandidateSpectrumLibFactory,
 * use both parent mass and peaks as filter, the reason the x
 * near the class name. Here rather than hashing peptides by mass
 * we take a protein centric approach that should be more easily 
 * scale to all peptides scenario, candidate generation is done through
 * a linear scan through the database
 * @author Jian Wang
 *
 */
public class CandidateSpectrumLibFactoryX {
	protected int minCharge = 1;
	protected int maxCharge = 5;
	protected String proteinFile;
	protected StringBuffer proteins; //we construct the db by concatenating all proteins as one long string separated by a special characters
	protected double parentMassTolerance = 3;    
	protected String[] prefix = {"b"};
	protected String[] suffix = {"y"};
	protected int topPeaks = 20;  //parameters for peaks filter
	protected int matchTopPeaks = 3;
	protected double fragmentMassTolerance = 0.5;
	public CandidateSpectrumLibFactoryX(String proteinfile){
		this.proteinFile = proteinfile;
		readProteinFromFile();
		System.out.println("done reading proteins. DB size: " + this.proteins.length());
	}
	//TODO: we neeed to also retain which peptide belong to which parent protein
	//This is needed for two reasons
	//1) To obtain protein information for ID purpose and
	//2) To remember which peptides come from decoy database if such procedure is used
	//this can be done by probalby keeping track of which regions in the big proteins
	//belong to which proteins
	private void readProteinFromFile(){
		try{
			BufferedReader reader = new BufferedReader(new FileReader(this.proteinFile));
			String line = reader.readLine();
			proteins = new StringBuffer();
			while(line != null){
				if(!line.startsWith(">")){
					proteins.append(line);
				}else{
					proteins.append("*");
				}
				line = reader.readLine();
			}
			proteins.append("*"); //added a end token to make string manipulation easy
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public List<Peptide> getCandidatePeptides(Spectrum querySpectrum){
		Map<Double, Peak> table = new HashMap<Double, Peak>();
		List<Peak> queryPeaks = querySpectrum.getTopPeaks(this.topPeaks);
		for(int i = 0; i < queryPeaks.size(); i++){
			//System.out.println("Putting " + Math.round(queryPeaks.get(i).getMass()) + " in table");
			table.put(new Double(Math.round(queryPeaks.get(i).getMass())), queryPeaks.get(i));
		}
		double parentMass = querySpectrum.parentMass*querySpectrum.charge 
			- querySpectrum.charge*Mass.PROTON_MASS - Mass.WATER;
		System.out.println("parent mass is: " + parentMass);
		int length = this.proteins.length();
		int startInd = 0, endInd = 0;  //we keep two pointers to generate substrings
		double currentMass = 0.0;
		Set<Peptide> candidates = new HashSet();
		double[] chargedMasses = new double[this.maxCharge - this.minCharge + 1];
		//double[] bIons = new double[this.maxCharge - this.minCharge + 1];
		int matchCount = 0;
		int matchedSiteCount = 0;
		while(startInd < length && endInd < length && startInd <= endInd){   //forever
			//System.out.println("current mass: " + (currentMass + Mass.WATER) + " current peptide: " + this.proteins.substring(startInd, endInd));				
			for(int c = this.minCharge; c <= querySpectrum.charge; c++){
				int index = c-minCharge;
				chargedMasses[index] = (currentMass) / c;
			}
			matchCount = 0;
			if(this.proteins.charAt(endInd) == 'K'){
				matchedSiteCount++;
			}
			
			if(Math.abs(chargedMasses[0] - parentMass) < parentMassTolerance && matchedSiteCount >= 0){
				//System.out.println("match count: " + matchCount);
				//System.out.println("matched mass " + currentMass);
				double tempMass = currentMass;
				if(checkMatchedPeaks(startInd, endInd, currentMass, querySpectrum, parentMass, table)){
					generateCandidates(startInd, endInd, candidates, querySpectrum.charge);
				}
			}
			
			if(chargedMasses[0] < parentMass){
				if(endInd < length){
					//System.out.println("adding mass for: " + endInd);
					currentMass += Mass.getAAMass(this.proteins.charAt(endInd));
				}
				endInd++;
				continue;
			}
			if(chargedMasses[0] > parentMass){
//				System.out.println("mass too large");
				currentMass -= Mass.getAAMass(this.proteins.charAt(startInd));
				if(this.proteins.charAt(startInd) == 'K'){
					matchedSiteCount--;
				}
				startInd++;
			}
		}
			//System.out.println("i: " + startInd + "\t" + "j: " + endInd);	
		List<Peptide> peps = new ArrayList();
		peps.addAll(candidates);
		return peps;
   }	
	
	private boolean checkMatchedPeaks(int startInd, int endInd, 
			double currentMass, Spectrum querySpectrum, 
			double parentMass, Map<Double, Peak> peakTable){
		double tempMass = currentMass;
		int matchCount = 0;
		for(int index = endInd-1; index > startInd; index--){				
			tempMass -= Mass.getAAMass(this.proteins.charAt(index));
			for(int charge = minCharge; charge <= querySpectrum.charge; charge++){
				double b = (tempMass + Mass.getIonMod("b") + (charge-1)*Mass.PROTON_MASS)/charge;
				double y = (parentMass - tempMass + Mass.getIonMod("y")+ (charge-1)*Mass.PROTON_MASS)/charge;
				//System.out.print(b + "| ");
				Double key = new Double(Math.round(b));
				Double key2 = new Double(Math.round(y));
				if(peakTable.containsKey(key) || peakTable.containsKey(key2)){  //let's try b-ions and y-ions first
					matchCount++;
				}
			}
			if(matchCount >= this.matchTopPeaks){
				return true;
			}					
		}
		return false;
	}
	
	protected void generateCandidates(int begin, int end, Set<Peptide> candidates, int charge){
		Peptide cand = new Peptide(this.proteins.substring(begin, end), charge);
		candidates.add(cand);
	}
	
	//TODO: Future implementation should include batch look up as well to
	//speed up the searches as well as amortized the time spent for computing
	//all parente mass and b/y frangments
	public List<Peptide> getCandidatePeptides(List<Spectrum> querySpectrum){
		return null;
	}
	
	public static boolean checkPassFilter(String peptide, List<Peptide> filtered){
		String p = peptide.split("\\.")[0];
		for(int i = 0; i < filtered.size(); i++){
			Peptide curr = filtered.get(i);
			if(p.equals(curr.getPeptide())){
				return true;
			}
		}
		return false;
	}
	

	
	public static void testCandidateSpectrumLibFactoryX(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		lib1.removeModSpectra();
		int testSetSize = 1000;
		System.out.println("Done loading library");
		List<Spectrum> testset = new ArrayList<Spectrum>();
		for(int i = 0; i < testSetSize; i++){
			Spectrum testSpectrum = lib1.getRandomSpectrum();
			testset.add(testSpectrum);
		}
		lib1 = null;
		String proteinFile = "..\\mixture_linked\\Ecoli_genome.fasta";
		String peptideFile = "..\\mixture_linked\\Ecoli_allpeptides.txt";
		CandidateSpectrumLibFactoryX xfactory = new CandidateSpectrumLibFactoryX(proteinFile);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < testSetSize; i++){
			Spectrum testSpectrum = testset.get(i);
			List<Peptide> candidates = xfactory.getCandidatePeptides(testSpectrum);
			System.out.print("Query " + testSpectrum.peptide +  " Parentmass filter found: " + candidates.size() + " candidates, include target: ");
			System.out.println(checkPassFilter(testSpectrum.peptide, candidates));
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void main(String[] args){
		testCandidateSpectrumLibFactoryX();
	}
	
}
