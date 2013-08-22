package org.Spectrums;

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;


	public class LinkedCandidateSpectrumLibX extends CandidateSpectrumLibFactoryX{
		private int minCommonCharge = 1;
		private int maxCommonCharge = 2;
		private int minLinkedCharge = 2;
		private int maxLinkedCharge = 6;
		
		public LinkedCandidateSpectrumLibX(String proteinFile){
			super(proteinFile);
		}
		
		public List<String> filterByPeaks(Spectrum querySpectrum){
//			System.out.println("starting peak filter");
			Map<Double, List<Peak>> table = new HashMap<Double, List<Peak>>();
			generateQueryPeaks(querySpectrum, table);
			double parentMass = querySpectrum.parentMass*querySpectrum.charge 
				- querySpectrum.charge*Mass.PROTON_MASS - Mass.WATER;
			System.out.println("parent mass is: " + parentMass);
			int length = this.proteins.length();
			int startInd = 0, endInd = startInd+1;  //we keep two pointers to generate substrings
			//Set<String> candidates = new HashSet();
			Map<Pair, Integer> candidates = new HashMap<Pair, Integer>();
			while(startInd < length-1 && startInd <= endInd){
				double currentMass = Mass.getAAMass(proteins.charAt(startInd));
				int matchCount = 0;
				Set<Peak> matchedPeaks = new HashSet();
				while(currentMass <  parentMass - 00 && endInd < length){
					for(int charge = 1; charge < 2; charge++){
						double b = (currentMass + Mass.getIonMod("b") + (charge-1)*Mass.PROTON_MASS)/charge;
						//System.out.println("b: " + b + " for: " + proteins.substring(startInd, endInd));
						if(table.containsKey(discretize(b))){
							Peak p = table.get(discretize(b)).get(0);
							double k = p.getMass();
							if(!matchedPeaks.contains(p)){
								matchedPeaks.add(p);								
								matchCount++;
								//System.out.println("match count is: " + matchCount);
							}
						}
						if(matchCount > 0){
							if(endInd - startInd >= 5 ){
								int i = startInd > 0 ? startInd-1 : 0;
								int j = endInd > 0 ? endInd-1 : 0;
								if((proteins.charAt(i) == 'K' || proteins.charAt(i) == 'R')
										|| (proteins.charAt(j) == 'K' || proteins.charAt(j) == 'R')){
									incrementCounter(candidates, startInd, endInd, matchCount);
									//System.out.println("adding count " + proteins.substring(startInd, endInd) + " with " + matchCount);
								}
							}
						}
					}
					currentMass += Mass.getAAMass(proteins.charAt(endInd++));
				}
				startInd++;
				endInd = startInd+1;
				
			}
			
			startInd = proteins.length()-1;
			endInd = startInd - 1;
			while(startInd > 0 && startInd >= endInd){
				double currentMass = Mass.getAAMass(proteins.charAt(startInd));
				int matchCount = 0;
				Set<Peak> matchedPeaks = new HashSet();
				while(currentMass <  parentMass - 00 && endInd >= 0){
					for(int charge = 1; charge < 2; charge++){
						double y = (currentMass + Mass.getIonMod("y") + (charge-1)*Mass.PROTON_MASS)/charge;
						//System.out.println("y: " + y + " for: " + proteins.substring(endInd+1, startInd+1));
						if(table.containsKey(discretize(y))){
							//System.out.println("matching");
							Peak p = table.get(discretize(y)).get(0);
							double k = p.getMass();
							if(!matchedPeaks.contains(p)){
								matchedPeaks.add(p);
								matchCount++;
							}
						}
						if(matchCount > 0){
							if(startInd - endInd >= 5 ){
								int i = startInd < length ? startInd : length-1;
								int j = endInd < length ? endInd : length-1;
								if((proteins.charAt(i) == 'K' || proteins.charAt(i) == 'R')
										|| (proteins.charAt(j) == 'K' || proteins.charAt(j) == 'R')){
									incrementCounter(candidates, endInd+1, startInd+1, matchCount);
									//System.out.println("adding count " + proteins.substring(endInd+1, startInd+1) + " with "  + matchCount);
								}
							}
						}
					}
					currentMass += Mass.getAAMass(proteins.charAt(endInd--));
				}
				startInd--;
				endInd = startInd-1;
			}
			
			List<String> peps = new ArrayList();
			getCandidates(candidates, peps);
			return peps;
	    }
		
		public List<String> filterByPeaks3(Spectrum querySpectrum){
			System.out.println("starting peak filter");
			Map<Double, List<Peak>> table = new HashMap<Double, List<Peak>>();
			generateQueryPeaks(querySpectrum, table);
			double parentMass = querySpectrum.parentMass*querySpectrum.charge 
				- querySpectrum.charge*Mass.PROTON_MASS - Mass.WATER;
			System.out.println("parent mass is: " + parentMass);
			int length = this.proteins.length();
			int startInd = 0, endInd = startInd+1;  //we keep two pointers to generate substrings
			Set<String> candidates = new HashSet();
			//Map<Pair, Integer> candidates = new HashMap<Pair, Integer>();
			while(startInd < length-1 && startInd <= endInd){
				double currentMass = Mass.getAAMass(proteins.charAt(startInd));
				int matchCount = 0;
				Set<Peak> matchedPeaks = new HashSet();
				while(currentMass <  parentMass - 00 && endInd < length){
					double b = (currentMass + Mass.getIonMod("b"));
					//System.out.println("b: " + b + " for: " + proteins.substring(startInd, endInd));
					if(table.containsKey(discretize(b))){
						Peak p = table.get(discretize(b)).get(0);
						double k = p.getMass();
						if(!matchedPeaks.contains(p)){
							matchedPeaks.add(p);								
							matchCount++;
							//System.out.println("match count is: " + matchCount);
						}
					}
					if(matchCount >= this.matchTopPeaks){
						System.out.println("peptide pass filter");
						if(endInd - startInd >= 5 ){
							int i = startInd > 0 ? startInd-1 : 0;
							int j = endInd > 0 ? endInd-1 : 0;
							if((proteins.charAt(i) == 'K' || proteins.charAt(i) == 'R')
									|| (proteins.charAt(j) == 'K' || proteins.charAt(j) == 'R')){
								//candidates.add(proteins.substring(startInd, endInd));
							}
						}
					}
					
					currentMass += Mass.getAAMass(proteins.charAt(endInd++));
				}
				//backward pass
				endInd -= 1;	
				currentMass = Mass.getAAMass(proteins.charAt(endInd));
				while(currentMass < parentMass - 00 && endInd > 0){
					double y = (currentMass + Mass.getIonMod("y"));
					//System.out.println("b: " + b + " for: " + proteins.substring(startInd, endInd));
					if(table.containsKey(discretize(y))){
						Peak p = table.get(discretize(y)).get(0);
						double k = p.getMass();
						if(!matchedPeaks.contains(p)){
							matchedPeaks.add(p);								
							matchCount++;
							//System.out.println("match count is: " + matchCount);
						}
						if(matchCount >= this.matchTopPeaks){
							if(endInd - startInd >= 5 ){
								int i = startInd > 0 ? startInd-1 : 0;
								int j = endInd > 0 ? endInd : 0;
								//if((proteins.charAt(i) == 'K' || proteins.charAt(i) == 'R')
								//		|| (proteins.charAt(j) == 'K' || proteins.charAt(j) == 'R')){
									System.out.println("adding peptides");
									candidates.add(proteins.substring(startInd, endInd+1));
									
								//}
							}
						}
					}
					currentMass += Mass.getAAMass(proteins.charAt(endInd--));
				}
				
				startInd++;
				endInd = startInd+1;	
			}
			List<String> peps = new ArrayList();
			peps.addAll(candidates);
			return peps;
		}
		
		//uncombined version
		public List<String> filterByPeaks2(Spectrum querySpectrum){
//			System.out.println("starting peak filter");
			Map<Double, List<Peak>> table = new HashMap<Double, List<Peak>>();
			generateQueryPeaks(querySpectrum, table);
			double parentMass = querySpectrum.parentMass*querySpectrum.charge 
				- querySpectrum.charge*Mass.PROTON_MASS - Mass.WATER;
			System.out.println("parent mass is: " + parentMass);
			int length = this.proteins.length();
			int startInd = 0, endInd = startInd+1;  //we keep two pointers to generate substrings
			Set<String> candidates = new HashSet();
			//Map<Pair, Integer> candidates = new HashMap<Pair, Integer>();
			while(startInd < length-1 && startInd <= endInd){
				double currentMass = Mass.getAAMass(proteins.charAt(startInd));
				int matchCount = 0;
				Set<Peak> matchedPeaks = new HashSet();
				while(currentMass <  parentMass - 00 && endInd < length){
					for(int charge = 1; charge < 2; charge++){
						double b = (currentMass + Mass.getIonMod("b") + (charge-1)*Mass.PROTON_MASS)/charge;
						//System.out.println("b: " + b + " for: " + proteins.substring(startInd, endInd));
						if(table.containsKey(discretize(b))){
							Peak p = table.get(discretize(b)).get(0);
							double k = p.getMass();
							if(!matchedPeaks.contains(p)){
								matchedPeaks.add(p);								
								matchCount++;
								//System.out.println("match count is: " + matchCount);
							}
						}
						if(matchCount >= this.matchTopPeaks){
							if(endInd - startInd >= 5 ){
								int i = startInd > 0 ? startInd-1 : 0;
								int j = endInd > 0 ? endInd-1 : 0;
								if((proteins.charAt(i) == 'K' || proteins.charAt(i) == 'R')
										|| (proteins.charAt(j) == 'K' || proteins.charAt(j) == 'R')){
									candidates.add(proteins.substring(startInd, endInd));
								}
							}
						}
					}
					currentMass += Mass.getAAMass(proteins.charAt(endInd++));
				}
				startInd++;
				endInd = startInd+1;
				
			}
			
			startInd = proteins.length()-1;
			endInd = startInd - 1;
			while(startInd > 0 && startInd >= endInd){
				double currentMass = Mass.getAAMass(proteins.charAt(startInd));
				int matchCount = 0;
				Set<Peak> matchedPeaks = new HashSet();
				while(currentMass <  parentMass - 00 && endInd >= 0){
					for(int charge = 1; charge < 2; charge++){
						double y = (currentMass + Mass.getIonMod("y") + (charge-1)*Mass.PROTON_MASS)/charge;
						//System.out.println("y: " + y + " for: " + proteins.substring(endInd+1, startInd+1));
						if(table.containsKey(discretize(y))){
							//System.out.println("matching");
							Peak p = table.get(discretize(y)).get(0);
							double k = p.getMass();
							if(!matchedPeaks.contains(p)){
								matchedPeaks.add(p);
								matchCount++;
							}
						}
						if(matchCount >= this.matchTopPeaks){
							if(startInd - endInd >= 5 ){
								int i = startInd < length ? startInd : length-1;
								int j = endInd < length ? endInd : length-1;
								if((proteins.charAt(i) == 'K' || proteins.charAt(i) == 'R')
										|| (proteins.charAt(j) == 'K' || proteins.charAt(j) == 'R')){
									candidates.add(this.proteins.substring(endInd+1, startInd+1));
								}
							}
						}
					}
					currentMass += Mass.getAAMass(proteins.charAt(endInd--));
				}
				startInd--;
				endInd = startInd-1;
			}
			
			List<String> peps = new ArrayList();
			peps.addAll(candidates);
			return peps;
	    }
		//ALKAWSVAR & VHKECCHGDLLECADDRADLAK
		private void generateQueryPeaks(Spectrum query, Map<Double, List<Peak>> table){
			List<Peak> queryPeaks = query.getTopPeaks(this.topPeaks);
			double parentMass = query.parentMass*query.charge - Mass.WATER;
			for(int i = 0; i < queryPeaks.size(); i++){
				//System.out.println("Putting " + Math.round(queryPeaks.get(i).getMass()) + " in table");
				Peak p = queryPeaks.get(i);
				for(int charge = 1; charge <= 2; charge++){
					Double key = new Double(Math.round(p.getMass()*charge - (charge-1)*Mass.PROTON_MASS));
					//System.out.println("putting key " + key);
					insertPeakTable(key, p, table);
				}
				for(int charge = 2; charge < 5; charge++){
					double complement = (parentMass - p.getMass()*charge) - Mass.getIonMod("b") + Mass.getIonMod("y");  //assume peak is a b ions
					//System.out.println("complement is: " + complement);
					insertPeakTable(new Double(Math.round(complement)), p, table);
					insertPeakTable(new Double(Math.round(complement-1.0)), p, table);
					insertPeakTable(new Double(Math.round(complement+1.0)), p, table);
					double complement2 = (parentMass - p.getMass()*charge) - Mass.getIonMod("y") + Mass.getIonMod("b"); //assume peak is a y ions 
					//System.out.println("complement2 is: " + complement);
					insertPeakTable(new Double(Math.round(complement2)), p, table);
					insertPeakTable(new Double(Math.round(complement2-1.0)), p, table);
					insertPeakTable(new Double(Math.round(complement2+1.0)), p, table);

				}
			}
		}
		
		
		private void insertPeakTable(Double key, Peak p, Map<Double, List<Peak>> table){
			List<Peak> pList;
			if(table.containsKey(key)){
				pList = table.get(key);
			}else{
				pList = new ArrayList<Peak>();
			}
			pList.add(p);
			table.put(key, pList);
		}
		
		private void getCandidates(Map<Pair, Integer> candidates, List<String> candidateHits){
			Set<String> uniqueCandidates = new HashSet<String>();
			for(Iterator<Pair> it = candidates.keySet().iterator(); it.hasNext();){
				Pair curr = it.next();
				if(candidates.get(curr) >= this.matchTopPeaks){
					uniqueCandidates.add(this.proteins.substring(((Integer)curr.getFirst()).intValue(), 
							((Integer)curr.getSecond()).intValue()));
				}
			}
			candidateHits.addAll(uniqueCandidates);
		}
		
		private void incrementCounter(Map<Pair, Integer> candidates, int startInd, int endInd, int increment){
			Pair key = new Pair(new Integer(startInd), new Integer(endInd));
			if(candidates.containsKey(key)){
				Integer count = candidates.get(key);
				count=count+increment;
				candidates.put(key, count);
				//System.out.println("incrememnt count");
				//System.out.println("after increment count is: " + candidates.get(key));
			}else{
				//System.out.println("Generating new entry");
				candidates.put(key, new Integer(increment));
				//Pair temp = new Pair(new Integer(startInd), new Integer(endInd));
				//System.out.println("Does table store entry successfully: " + candidates.containsKey(temp));
			}
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
		
		public static int checkPassFilter(String peptide1, String peptide2, List<String> filtered){
			int matchCount = 0;
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
		

		public static boolean checkPassFilter(String peptide1, List<String> filtered){
			int matchCount = 0;
			peptide1 = peptide1.split("\\.")[0];
			for(int i = 0; i < filtered.size(); i++){
				String curr = filtered.get(i);
				if(curr.equals(peptide1)){
					matchCount++;
				}
			}
			return matchCount > 0;
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
			LinkedCandidateSpectrumLibX xfactory = new LinkedCandidateSpectrumLibX(proteinFile);
			long start = (new GregorianCalendar()).getTimeInMillis();
			for(int i = 0; i < testSetSize; i++){
				Spectrum testSpectrum = testset.get(i);
				xfactory.matchTopPeaks = 2;
				xfactory.topPeaks = 40;
				List<String> candidates = xfactory.filterByPeaks(testSpectrum);
				CrossLinker linker = new CrossLinker(Mass.DSPLINKER_MASS, new char[]{'K'});
				LinkedCandidateSpectrumLibFactory f = new LinkedCandidateSpectrumLibFactory(candidates, linker);
				List<Peptide> candidates2 = f.getLinkedCandidateByMass(testSpectrum.parentMass+600, 0.2, testSpectrum.charge);
				System.out.print("Query " + testSpectrum.peptide +  " ions mass filter found: " + candidates.size() + " candidates, include target: ");
				System.out.println(checkPassFilter(testSpectrum.peptide, candidates));
				System.out.println("Parent mass filter found: " + candidates2.size() + " candidates");
			}
			System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

		}
		public static void testLinkCandidateSpectrumLibX(){
//			String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
//			SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
//			lib1.removeModSpectra();
			String spectrumFile = "..\\mixture_linked\\spectrums_raw.mgf";
			SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
			System.out.println("Done loading library");
			List<Spectrum> testset = new ArrayList<Spectrum>();
			testset.addAll(lib1.getAllSpectrums());
			lib1 = null;
			//String proteinFile = "..\\mixture_linked\\Ecoli_genome.fasta";
			String proteinFile = "..\\mixture_linked\\protein_linkedpeptides.fasta";
			String peptideFile = "..\\mixture_linked\\Ecoli_allpeptides.txt";
			LinkedCandidateSpectrumLibX xfactory = new LinkedCandidateSpectrumLibX(proteinFile);
			long start = (new GregorianCalendar()).getTimeInMillis();
			for(int i = 0; i < 81; i++){
				Spectrum testSpectrum = testset.get(i);
				List<String> candidates = xfactory.filterByPeaks3(testSpectrum);
				xfactory.matchTopPeaks = 2;
				CrossLinker linker = new CrossLinker(Mass.DSPLINKER_MASS, new char[]{'K'});
				LinkedCandidateSpectrumLibFactory fact = new LinkedCandidateSpectrumLibFactory(candidates, linker);
				fact.setMaxCharge(1);
				System.out.println("Query " + testSpectrum.peptide +  " b/y fragment mass filter found: " + candidates.size() + " candidates, include target: ");
				List<Peptide> filtered = fact.getLinkedCandidateByMass(testSpectrum.parentMass, 1.0, testSpectrum.charge);
				System.out.print("Query " + testSpectrum.peptide +  " Parentmass filter found: " + filtered.size() + " candidates, include target: ");
				String[] peptides = testSpectrum.peptide.split(" & ");
				System.out.println(checkPassFilter(peptides[0], peptides[1].substring(0, peptides[1].length()-2), candidates));
//				System.out.println("")
			}
			System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

		}
		
		public static void main(String[] args){
			//testCandidateSpectrumLibFactoryX();
			testLinkCandidateSpectrumLibX();
		}
}
