package org.Spectrums;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections.MultiMap;
import org.apache.commons.collections.map.MultiValueMap;

public class DecoySpectrumGenerator {
	public static char[] StaticResidues = {'K', 'R', 'P'}; 
	private static String DECOY="DECOY";
	private Spectrum s;
	public double tolerance =0.05;
	public DecoySpectrumGenerator(){
	}
	
	public Spectrum generateDecoy1(Spectrum s){
		this.s = s;
		Spectrum d = new Spectrum();
		d.spectrumName = s.spectrumName;
		d.scanNumber = s.scanNumber;
		d.charge = this.s.charge;
		d.parentMass = this.s.parentMass;
		Peptide pep = new Peptide(s.peptide, s.charge);
		d.peptide = shufflePeptide(pep.getPeptide());
		Peptide decoyPep = new Peptide(d.peptide, s.charge);
		TheoreticalSpectrum.prefixIons = new String[]{"b", "b-H20", "b-NH3"};//, "b(iso)"};
		TheoreticalSpectrum.suffixIons = new String[]{"y", "y-H20", "y-NH3"};//, "y(iso)"};
		TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
		TheoreticalSpectrum.addIsotopicPeaks(t, 1);
		SimpleMatchingGraph g = t.getMatchGraph(s, tolerance);
		MultiMap map = generateAnnotationMap(g);
		TheoreticalSpectrum t2 = new TheoreticalSpectrum(decoyPep);
		TheoreticalSpectrum.addIsotopicPeaks(t2, 1);
		List<Peak> pList = new ArrayList<Peak>();
		for(Iterator<Peak> it = t2.getPeak().iterator(); it.hasNext();){
			LabelledPeak lp = (LabelledPeak)it.next();
			String key = lp.getPeakType();
			if(map.containsKey(key)){
				for(Iterator<Peak> it2 = ((Collection)map.get(key)).iterator(); it2.hasNext();){
					Peak oldPeak = it2.next();
					Peak newPeak = new Peak(lp.getMass(), oldPeak.getIntensity());
					pList.add(newPeak);
				}
			}
		}
		//System.out.println(s.spectrumName + "\t" + s.peptide + " number of annotated: " + pList.size());
		addUnannotatedPeaks(g, pList);
		//s.getPeak().removeAll(getUnAnnotatedPeak(g));
		Collections.sort(pList, PeakMassComparator.comparator);
		d.setPeaks(pList);
		d.spectrumName = "DECOY_" + d.spectrumName;
		d.protein = "DECOY_" + s.protein;
		//d.peptide = "X_"+shufflePeptide(s.peptide);
		return d;
	}
	
	public Spectrum generateDecoy(Spectrum s){
		this.s = s;
		Spectrum d = new Spectrum();
		d.spectrumName = s.spectrumName;
		d.scanNumber = s.scanNumber;
		d.charge = this.s.charge;
		d.parentMass = this.s.parentMass;
		d.rt = s.rt;
		Peptide pep = new Peptide(s.peptide, s.charge);
		d.peptide = shufflePeptide(pep.getPeptide());
		Peptide decoyPep = new Peptide(d.peptide, s.charge);
		TheoreticalSpectrum.prefixIons = new String[]{"b", "b-H20", "b-NH3"};//, "b(iso)"};
		TheoreticalSpectrum.suffixIons = new String[]{"y", "y-H20", "y-NH3"};//, "y(iso)"};
		TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
		TheoreticalSpectrum.addIsotopicPeaks(t, 1);
		SimpleMatchingGraph g = t.getMatchGraph(s, tolerance);
		MultiMap map = generateAnnotationMap(g);
		TheoreticalSpectrum t2 = new TheoreticalSpectrum(decoyPep);
		TheoreticalSpectrum.addIsotopicPeaks(t2, 1);
		List<Peak> pList = new ArrayList<Peak>();
		Map<String, LabelledPeak> map2 = new HashMap<String, LabelledPeak>();
		for(Iterator<Peak> it = t2.getPeak().iterator(); it.hasNext();){
			LabelledPeak lp = (LabelledPeak)it.next();
			String key = lp.getPeakType();
			map2.put(key, lp);
		}
		
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			LabelledPeak closestP = null;
			double smallestErr = 10000;
			for(int i = 0; i < neighs.size(); i++){
				LabelledPeak currPeak = (LabelledPeak)neighs.get(i);
				double currDiff = Math.abs(currPeak.getMass()-p.getMass());
				closestP = currDiff < smallestErr ? currPeak : closestP;
				smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
			}
			if(closestP != null){
				LabelledPeak lp = closestP;
				String key = lp.getPeakType();
				if(map2.containsKey(key)){
					LabelledPeak newlp = map2.get(key);
					Peak newPeak = new Peak(newlp.getMass(), p.getIntensity());
					pList.add(newPeak);
					//System.out.println("DECOY Peaks: " + newPeak + "\t" + newlp);
				}
			}
			
		}
		
		
		//System.out.println(s.spectrumName + "\t" + s.peptide + " number of annotated: " + pList.size());
		addUnannotatedPeaks(g, pList);
		//s.getPeak().removeAll(getUnAnnotatedPeak(g));
		Collections.sort(pList, PeakMassComparator.comparator);
		d.setPeaks(pList);
		d.spectrumName = "DECOY_" + d.spectrumName;
		d.protein = "DECOY_" + s.protein;
		//d.peptide = "X_"+shufflePeptide(s.peptide);
		return d;
	}
	
	
	
	//generate a spectrum for target, but correct all observed masses
	//to theoretical masses position
	
	public Spectrum generateTarget(Spectrum s){
		return generateTarget(s, new String[]{"b", "b-H20", "b-NH3"}, new String[]{"y", "y-H20", "y-NH3"});
	}
	
	public Spectrum generateTarget(Spectrum s, String[] prefixIons, String[] suffixIons){
		this.s = s;
		Spectrum t = new Spectrum();
		t.spectrumName = s.spectrumName;
		t.scanNumber = s.scanNumber;
		t.charge = this.s.charge;
		t.parentMass = this.s.parentMass;
		t.peptide = s.peptide;
		t.protein = s.protein;
		t.rt = s.rt;
		Peptide pep = new Peptide(s.peptide, s.charge); 
		//TheoreticalSpectrum.prefixIons = prefixIons;
		//TheoreticalSpectrum.suffixIons = suffixIons;
		TheoreticalSpectrum th = new TheoreticalSpectrum(pep, prefixIons, suffixIons);
		TheoreticalSpectrum.addIsotopicPeaks(th, 1);
		SimpleMatchingGraph g = th.getMatchGraph(s, tolerance);
		List<Peak> pList = new ArrayList<Peak>();
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			LabelledPeak closestP = null;
			double smallestErr = 10000;
			for(int i = 0; i < neighs.size(); i++){
				LabelledPeak currPeak = (LabelledPeak)neighs.get(i);
				double currDiff = Math.abs(currPeak.getMass()-p.getMass());
				closestP = currDiff < smallestErr ? currPeak : closestP;
				smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
			}
			if(neighs.size() > 0){
				//System.out.println("annotated peaks\t" + p);
				//pList.add(p);
			}
			if(neighs.size() == 0){
				//System.out.println("unannotated peaks\t" + p);
			}
			if(closestP != null){
				Peak theo = new Peak(closestP.getMass(), p.getIntensity());
				//System.out.println("Target Peak: " + theo + "labelledP\t" + closestP);
				pList.add(theo);
			}
		}
		//System.out.println(s.spectrumName + "\t" + s.peptide + " number of annotated: " + pList.size() + "\t" + s.getPeak().size());
		addUnannotatedPeaks(g, pList);
		//s.getPeak().removeAll(getUnAnnotatedPeak(g));
		Collections.sort(pList, PeakMassComparator.comparator);
		t.setPeaks(pList);
		return t;
	}
	
	
	
	private void addUnannotatedPeaks(SimpleMatchingGraph g, List<Peak> pList){
		int annotated = pList.size();
		List<Peak> unAnnotated = new ArrayList();
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				//Peak newPeak = new Peak(p.getMass(), p.getIntensity());
				//Peak newPeak = new Peak(p.getMass(), 10);
				unAnnotated.add(p);
			}
		}
		//Collections.sort(unAnnotated, PeakIntensityComparator.comparator);
		//System.out.println("number of unannotated: " + (pList.size() - annotated) + " total: " + pList.size());
	}
	
	
	private List<Peak> getUnAnnotatedPeak(SimpleMatchingGraph g){
		List<Peak> unAnnotated = new ArrayList<Peak>();
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				unAnnotated.add(p);
			}
		}
		return unAnnotated;
	}
	
	private MultiMap generateAnnotationMap(TheoreticalSpectrum t, Spectrum s){
		return generateAnnotationMap(t.getMatchGraph(s, 0.5));
	}

	private MultiMap generateAnnotationMap(SimpleMatchingGraph mg){
		MultiMap annoMap = new MultiValueMap();
		for(Iterator<LabelledPeak> it = mg.vertexSet(mg.Theoretical).iterator(); it.hasNext();){
			LabelledPeak lp = it.next();
			List neighs = mg.getNeighbors(lp);
			for(Iterator<Peak> it2 = neighs.iterator(); it2.hasNext();){
				annoMap.put(lp.getPeakType(), it2.next());
			}
		}
		return annoMap;
	}
	//shuffle peptide
	public String shufflePeptide(String pep){
		String reverse = shuffle(pep);
		return reverse;
	}
	
	//shuffle string
	public String shuffle(String pep){
		StringBuffer shuffle = new StringBuffer(pep);
		List<Integer> index = new ArrayList<Integer>();  //list of index to shuffle
		//keep static position out of shuffle list
		for(int i = 0; i < pep.length(); i++){
			if(!isStaticResidue(pep.charAt(i))){
				index.add(i);
			}
		}
		if(index.size() == 0){
			return pep;
		}
		for(int i = 0; i < 50; i++){
			int j = (int)Math.floor(Math.random()*index.size());
			int k = (int)Math.floor(Math.random()*index.size());
			j = index.get(j);
			k = index.get(k);
			char c = shuffle.charAt(j);
			shuffle.setCharAt(j, shuffle.charAt(k));
			shuffle.setCharAt(k, c);
		}
		return shuffle.toString();
	}
	
//	public String shuffle(String pep){
//		StringBuffer shuffle = new StringBuffer(pep);
//		//keep static position out of shuffle list
//		for(int i = 0; i < 50; i++){
//			int j = (int)Math.floor(Math.random()*(pep.length()-2));
//			int k = (int)Math.floor(Math.random()*(pep.length()-2));
//			char c = shuffle.charAt(j);
//			shuffle.setCharAt(j, shuffle.charAt(k));
//			shuffle.setCharAt(k, c);
//		}
//		return shuffle.toString();
//	}

	
	public boolean isStaticResidue(char c){
		for(int i = 0; i < this.StaticResidues.length; i++){
			if(c == this.StaticResidues[i]){
				return true;
			}
		}
		return false;
	}
	
	public static void testShufflePeptide(){
		DecoySpectrumGenerator d = new DecoySpectrumGenerator();
		String peptide = "LGRQQADLLSVEDR";
		System.out.println("peptide is: " + peptide + "\tshuffle is: " + d.shuffle(peptide));
	}
	
	public static void testShufflePeptides(){
		String filename = "..\\mixture_linked\\database\\rabbit_proteasome_uniprot_rabbitseqs_msgfdbIDs.txt";
		DecoySpectrumGenerator d = new DecoySpectrumGenerator();
		List<String> peptides = Utils.FileIOUtils.createListFromFile(filename);
		int noDecoy=0;
		for(Iterator<String> it = peptides.iterator(); it.hasNext();){
			String pep = it.next();
			//System.out.println("pep is: " + pep);
			String decoy = d.shuffle(pep);
			if(decoy != null){	
				System.out.println(d.shuffle(pep));
			}else{
				noDecoy++;
			}
		}
		//System.out.println("Cannot generate decoy: " + noDecoy);
	}
	
	public static void testShuffleSpectra(String spectrumFile, String decoyFileName, double tolerance){
		//spectrumFile = "..\\mixture_linked\\SWATH-MSPLITv1.0/Testing/test_buildlib/Spectral_library.txt";
		//tolerance = 0.33;
		LargeSpectrumLibIterator<Spectrum> iter = new LargeSpectrumLibIterator(spectrumFile);
		int index = spectrumFile.lastIndexOf('.');
		if(decoyFileName == null || decoyFileName.length() < 1){
			decoyFileName = spectrumFile.substring(0, index) + "_plusDecoy.mgf";
		}
		try{
			BufferedWriter bo = new BufferedWriter(new FileWriter(decoyFileName));
			int fail=0;
			for(;iter.hasNext();){
				Spectrum s = iter.next();
				//System.out.println("Spectrum is : " + s.spectrumName);
				//s.shiftSpectrumPPM(100);
				//s.windowFilterPeaks(10, 25);
				//if(s.peptide.contains("+") || s.peptide.contains("(") || s.peptide.contains("-")){
					//continue;
				//}
				if(s.peptide.contains("+")){
					//continue;
				}
				DecoySpectrumGenerator generator = new DecoySpectrumGenerator();
				generator.tolerance = tolerance;
				//double signals = SpectrumUtil.getSNR(s, tolerance);
				//s.filterPeaksByRankScore(120);
				//if(signals < 8){
					//continue;
				//}
				Spectrum target = generator.generateTarget(s);
				if(target.getPeak().size() < 10){
					continue;
				}
				
				Spectrum decoy = generator.generateDecoy1(s);
				
				double sim = s.cosine(decoy, tolerance);
				int count=0;
				while(sim > 0.65){
					decoy = generator.generateDecoy(s);
					sim = s.cosine(decoy, tolerance);
					count++;
					if(count > 50){
						fail++;
						System.out.println("cannot generated decoy " + s.spectrumName + "\t" + s.peptide);
						decoy = null;
						break;
					}
				}
				//only keep significant peaks in library
				//s.shiftSpectrumPPM(-100);
				//bo.write(s.toString());
				if(decoy != null){
					//decoy.shiftSpectrumPPM(-100);
					bo.write(target.toString());
					//bo.write(s.toString());
					bo.write("\n");
					//System.out.println("target similarity to theo-target: " + s.cosine(target, tolerance));
					//System.out.println("decoy similarity to target: " + s.cosine(decoy, tolerance));
					//System.out.println("decoy projected-similarity to target: " + s.projectedCosine(decoy, 0.05));
					bo.write(decoy.toString());
					bo.write("\n");
				}
			}
			System.out.println("cannot generated decoy for " + fail + " spectra in the library");
			bo.flush();
			bo.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		}
		
	}
	
	public static void testDecoyGenerator(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		//String spectrumFile2 = "..\\mixture_linked\\yeast_annotated_spectra.mgf";
		String spectrumFile2 = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest_080109235858.mgf";
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MSP");
		SpectrumLib lib2 = new SpectrumLib(spectrumFile2, "MGF");
		//lib2 = lib2.Divide();
		lib.removeModSpectra();
		List<Spectrum> specList = lib.getAllSpectrums();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			DecoySpectrumGenerator generator = new DecoySpectrumGenerator();
			Spectrum decoy = generator.generateDecoy(s);
			//System.out.println(s);
			//System.out.println(decoy);
			lib.addSpectrum(decoy);
		}
		lib.scaleSpectrumMass(0.9995);
		lib2.scaleSpectrumMass(0.9995);
		System.out.println("Done generating decoys");
		System.out.println("searching " + lib2.getAllSpectrums().size() + " spectra");
		List<Spectrum> specList2 = lib2.getAllSpectrums();
		long start = (new GregorianCalendar()).getTimeInMillis();
		SpectrumLib.runSearch(lib, specList2, 0.5, "abc", 3);
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void runSearchWithDecoyLib(String libFile, String queryFile){
		SpectrumLib lib = new SpectrumLib(libFile, "MSP");
		SpectrumLib lib2 = new SpectrumLib(queryFile, "MGF");
		//lib2 = lib2.Divide();
		lib.removeModSpectra();
		List<Spectrum> specList = lib.getAllSpectrums();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			DecoySpectrumGenerator generator = new DecoySpectrumGenerator();
			Spectrum decoy = generator.generateDecoy(s);
			//System.out.println(s);
			//System.out.println(decoy);
			lib.addSpectrum(decoy);
		}
		lib.scaleSpectrumMass(0.9995);
		lib2.scaleSpectrumMass(0.9995);
		System.out.println("Done generating decoys");
		System.out.println("searching " + lib2.getAllSpectrums().size() + " spectra");
		List<Spectrum> specList2 = lib2.getAllSpectrums();		
		SpectrumLib.runSearch(lib, specList2, 0.5, "abc", 3);
	}
	
	
	public static void main(String[] args){
		//testShufflePeptide();
		//testDecoyGenerator();
		testShuffleSpectra(args[0], args[1], Double.parseDouble(args[2]));
		//testShufflePeptides();
		//runSearchWithDecoyLib(args[0], args[1]);
	}
}
