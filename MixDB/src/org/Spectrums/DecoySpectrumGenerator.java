package org.Spectrums;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.collections.MultiMap;
import org.apache.commons.collections.map.MultiValueMap;

public class DecoySpectrumGenerator {
	public static char[] StaticResidues = {'K', 'R', 'P'}; 
	private static String DECOY="DECOY";
	private Spectrum s;
	public DecoySpectrumGenerator(){
	}
	
	public Spectrum generateDecoy(Spectrum s){
		this.s = s;
		Spectrum d = new Spectrum();
		d.spectrumName = s.spectrumName;
		d.scanNumber = s.scanNumber;
		d.charge = this.s.charge;
		d.parentMass = this.s.parentMass;
		d.peptide = shufflePeptide(s.peptide);
		Peptide pep = new Peptide(s.peptide, s.charge);
		Peptide decoyPep = new Peptide(d.peptide, s.charge);	
		TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
		SimpleMatchingGraph g = t.getMatchGraph(s, 0.5);
		MultiMap map = generateAnnotationMap(g);
		TheoreticalSpectrum t2 = new TheoreticalSpectrum(decoyPep);
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
		addUnannotatedPeaks(g, pList);
		Collections.sort(pList, PeakMassComparator.comparator);
		d.setPeaks(pList);
		d.spectrumName = "DECOY_" + d.spectrumName;
		//d.peptide = "X_"+shufflePeptide(s.peptide);
		return d;
	}
	
	private void addUnannotatedPeaks(SimpleMatchingGraph g, List<Peak> pList){
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				pList.add(p);
			}
		}
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
		String filename = "..\\mixture_linked\\database\\lib_ucsd3_plus_subpeptides.txt";
		DecoySpectrumGenerator d = new DecoySpectrumGenerator();
		List<String> peptides = Utils.FileIOUtils.createListFromFile(filename);
		for(Iterator<String> it = peptides.iterator(); it.hasNext();){
			String pep = it.next();
			//System.out.println(pep);
			System.out.println(d.shuffle(pep));
		}
	}
	
	public static void testShuffleSpectra(String spectrumFile){
		spectrumFile = "..\\mixture_linked\\NIST_human_IT_2010-01-14.mgf";
		LargeSpectrumLibIterator<Spectrum> iter = new LargeSpectrumLibIterator(spectrumFile);
		int index = spectrumFile.lastIndexOf('.');
		String decoyName = spectrumFile.substring(0, index) + "_plusDecoy.mgf";
		try{
			BufferedWriter bo = new BufferedWriter(new FileWriter(decoyName));
			for(;iter.hasNext();){
				Spectrum s = iter.next();
				if(s.peptide.contains("(")){
					continue;
				}
				DecoySpectrumGenerator generator = new DecoySpectrumGenerator();
				Spectrum decoy = generator.generateDecoy(s);
				bo.write(s.toString());
				bo.write("\n");
				bo.write(decoy.toString());
				bo.write("\n");
			}	
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
		testShuffleSpectra(args[0]);
		//testShufflePeptides();
		//runSearchWithDecoyLib(args[0], args[1]);
	}
}
