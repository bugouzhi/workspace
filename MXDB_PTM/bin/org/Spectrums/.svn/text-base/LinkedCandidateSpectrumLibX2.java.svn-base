package org.Spectrums;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class LinkedCandidateSpectrumLibX2 extends CandidateSpectrumLibFactoryX{
	private int minCommonCharge = 1;
	private int maxCommonCharge = 2;
	private int minLinkedCharge = 2;
	private int maxLinkedCharge = 6;
	
	public LinkedCandidateSpectrumLibX2(String proteinFile){
		super(proteinFile);
	}
	
	public List<String> filterByPeaks(Spectrum querySpectrum){
		Map<Double, Peak> table = new HashMap<Double, Peak>();
		List<Peak> queryPeaks = querySpectrum.getTopPeaks(this.topPeaks);
		Double maxMass = 0.0, minMass = 1000000.0;
		for(int i = 0; i < queryPeaks.size(); i++){
			//System.out.println("Putting " + Math.round(queryPeaks.get(i).getMass()) + " in table");
			Double key = new Double(Math.round(queryPeaks.get(i).getMass()));
			table.put(key, queryPeaks.get(i));
			maxMass = key > maxMass ? key : maxMass;
			minMass = key < minMass ? key : minMass;
		}
		double parentMass = querySpectrum.parentMass*querySpectrum.charge 
			- querySpectrum.charge*Mass.PROTON_MASS - Mass.WATER;
		System.out.println("parent mass is: " + parentMass);
		System.out.println("peak mass interval : " + minMass + " : " + maxMass);
		int length = this.proteins.length();
		int startInd = 0, endInd = 0;  //we keep two pointers to generate substrings
		Set<String> candidates = new HashSet();
		while(startInd < length && endInd < length){
			double currentMass = 0.0;
			int matchCount = 0;
			int commonMatch = 0;
			int linkedMatch = 0;
			int siteMatchCount = 0;
			while(endInd < length && currentMass <= (parentMass - 500)){  //rough estimate
				char aa = this.proteins.charAt(endInd);
				if(aa == 'K'){
					siteMatchCount++;
				}
				currentMass += Mass.getAAMass(aa);
				endInd++;
				//System.out.println("current mass is: " + currentMass + " : " + this.proteins.substring(startInd, endInd+1));
				if(currentMass >= 10 && currentMass/querySpectrum.charge <= (maxMass + 2)){ //only check table if it is within range
					//if current peptide is a common peaks
					if(siteMatchCount == 0){
						for(int charge = getMinCharge(querySpectrum.charge); charge <= getMaxCharge(querySpectrum.charge); charge++){
							double b = (currentMass + Mass.getIonMod("b") + (charge-1)*Mass.PROTON_MASS)/charge;
							System.out.println("currentmass: " + b + " for b: " + this.proteins.substring(startInd, endInd));
							double y = 0.0; //y-ions is linked so we do not compute them for now
							if(table.containsKey(discretize(b)) || table.containsKey(discretize(y))){
								System.out.println("matching");
								matchCount++;
							}
						}
					}else{ //now we gets to the tricky business of linked peptide 
						linkedMatch = 0;
						double tempMass = currentMass;
						double linkerOffset = parentMass - currentMass;
						for(int index = endInd-1; this.proteins.charAt(index) != 'K'; index--){								
							for(int charge = getLinkedMinCharge(querySpectrum.charge); charge <= getLinkedMaxCharge(querySpectrum.charge); charge++){
								double b = (tempMass + linkerOffset + Mass.getIonMod("b") + (charge-1)*Mass.PROTON_MASS)/charge;
								//System.out.println("linkedmass: " + b + " for b: " + this.proteins.substring(startInd, index+1));
								double y = 0.0;
								if(table.containsKey(discretize(b)) || table.containsKey(discretize(y))){
									linkedMatch++;
								}
							}
							tempMass -= Mass.getAAMass(proteins.charAt(index));
						}
						if((matchCount+linkedMatch) >= this.matchTopPeaks && (currentMass > 300)){//a simple heuristics to make sure the candidate peptide is not too short						 																	       
							candidates.add(this.proteins.substring(startInd, endInd));
						}
					}
				}
//				double diff = currentMass - parentMass;
//				if(matchCount >= this.matchTopPeaks && (currentMass > 300 && diff < 300) && siteMatchCount > 0){  //a simple heuristics to make sure the candidate peptide					 																	     //is not too short or too long   
//					candidates.add(this.proteins.substring(startInd, endInd));
//				}

			}
			//System.out.println("starting ind is: " + startInd);
			startInd++;
			endInd = startInd;
		}
		System.out.println("we have unique candidates: " + candidates.size());
		List<String> peps = new ArrayList();
		peps.addAll(candidates);
		return peps;
   }

	
	private Double discretize(double value){
		return new Double(Math.round(value));
	}
	
	private static int getMinCharge(int spectrumCharge){
		if(spectrumCharge < 4){
			return 1;
		}else{
			return 1;
		}
	}
	
	private static int getMaxCharge(int spectrumCharge){
		if(spectrumCharge == 1){
			return 1;
		}else if(spectrumCharge <= 4){
			return 2;
		}else{
			return 3;
		}
	}
	
	private static int getLinkedMinCharge(int spectrumCharge){
		if(spectrumCharge == 3){
			return 2;
		}else{
			return 3;
		}
	}
	
	private static int getLinkedMaxCharge(int spectrumCharge){
		if(spectrumCharge == 3){
			return 3;
		}else{
			return 4;
		}
	}
	
	public static boolean checkPassFilter(String peptide1, String peptide2, List<String> filtered){
		int matchCount = 0;
		for(int i = 0; i < filtered.size(); i++){
			String curr = filtered.get(i);
			if(curr.equals(peptide1)){
				matchCount++;
			}
			if(curr.equals(peptide2)){
				matchCount++;
			}
		}
		return matchCount > 0;
	}
	
	public static void testLinkCandidateSpectrumLibX(){
//		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
//		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
//		lib1.removeModSpectra();
		int testSetSize = 20;
		String spectrumFile = "..\\mixture_linked\\spectrums_raw.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
		System.out.println("Done loading library");
		List<Spectrum> testset = new ArrayList<Spectrum>();
		for(int i = 0; i < testSetSize; i++){
			Spectrum testSpectrum = lib1.getRandomSpectrum();
			testset.add(testSpectrum);
		}
		lib1 = null;
		//String proteinFile = "..\\mixture_linked\\Ecoli_genome.fasta";
		String proteinFile = "..\\mixture_linked\\ protein_linkedpeptides.fasta";
		String peptideFile = "..\\mixture_linked\\Ecoli_allpeptides.txt";
		LinkedCandidateSpectrumLibX2 xfactory = new LinkedCandidateSpectrumLibX2(proteinFile);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < testSetSize; i++){
			Spectrum testSpectrum = testset.get(i);
			List<String> candidates = xfactory.filterByPeaks(testSpectrum);
			xfactory.matchTopPeaks = 0;
			CrossLinker linker = new CrossLinker(Mass.DSPLINKER_MASS, new char[]{'K'});
			LinkedCandidateSpectrumLibFactory fact = new LinkedCandidateSpectrumLibFactory(candidates, linker);
			fact.setMaxCharge(1);
			System.out.println("Query " + testSpectrum.peptide +  " b/y fragment mass filter found: " + candidates.size() + " candidates, include target: ");
			List<Peptide> filtered = fact.getLinkedCandidateByMass(testSpectrum.parentMass, 1.0, testSpectrum.charge);
			System.out.print("Query " + testSpectrum.peptide +  " Parentmass filter found: " + filtered.size() + " candidates, include target: ");
			String[] peptides = testSpectrum.peptide.split(" & ");
			System.out.println(checkPassFilter(peptides[0], peptides[1].substring(0, peptides[1].length()-2), candidates));
//			System.out.println("")
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void main(String[] args){
		testLinkCandidateSpectrumLibX();
	}
}
