package Spectrum;
import java.util.ArrayList;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.io.*;
import org.apache.commons.collections.MultiMap;
import org.apache.commons.collections.map.MultiValueMap;

import IO.SortedSpectrumReader;
import Utils.LargeHashMap;


/**
 * Index spectrum library by precursor mass and store them in file
 * only load part of the library into memory as a cache. The size of cache
 * can be controlled by the number of object in LargeHashMap. Candidate spectra
 * from the library can be retrieved by their parentmass
 * @author Jian Wang
 *
 */

public class CandidateSpectrumGenerator {
	private String libraryFileName;
	private LargeHashMap table;
	private double avgTolerance=3.0;
	private double fragmentTolerance = 0.5;
	public CandidateSpectrumGenerator(){
	}
	
	public CandidateSpectrumGenerator(String libraryFile){
		this.libraryFileName = libraryFile;		
	}
	
	public void loadLibraryFromFile(String libraryObjectFile) {
		this.table = new LargeHashMap(libraryObjectFile);
		this.table.loadLibraryFromFile(libraryObjectFile);
	}
	
	/**
	 * Generated a indexed spectral library using file-base hashmap
	 * @param outfile
	 */
	public void generateIndexedLibrary(String outfile){
		SpectrumLib lib = new SpectrumLib(this.libraryFileName);
		lib.scaleSpectrumMass(0.9995);
		lib.toNormVector(fragmentTolerance*2, fragmentTolerance,  2000);
		lib.normIntensity();
		lib.toNormVector(fragmentTolerance*2, fragmentTolerance,  2000);
		Map specTable = indexSpecraByMass((List<Spectrum>)lib.getAllSpectrums());
		this.table = new LargeHashMap(outfile);
		this.table.buildTable(specTable);
		this.loadLibraryFromFile(outfile);
	}
	
	/**
	 * Generated a indexed spectral library using a stream of library spectra
	 * that are sorted by precursor masses, this method do not need to load
	 * all the spectra into memory. Instead it only store the spectra from the same
	 * table bin in memory, thus allow indexing of spectral library that are relatively
	 * large in memory efficient ways
	 * @param outfile
	 */
	public void generateIndexedLibraryFromSorted(String outfile){
		this.table = new LargeHashMap(outfile);
		SortedSpectrumReader iter = new SortedSpectrumReader(this.libraryFileName);
		int prevKey = 0;
		List<Spectrum> specList = new ArrayList();
		System.out.println("Starting indexing");
		int count=0;
		while(iter.hasNext()){
			Spectrum s = (Spectrum)iter.next();
			s.toNormVector(fragmentTolerance*2,  fragmentTolerance, 2000);
			s.sqrtSpectrum();
			s.toNormVector(fragmentTolerance*2,  fragmentTolerance, 2000);
			int currentKey = getKey(s);
			if(prevKey == 0 || currentKey == prevKey){  
				specList.add(s);
			}else{
				this.table.put(currentKey, specList);
				specList.clear();
			}
			prevKey = currentKey;			
			count++;
			if(count % 1000 == 0){
				System.out.println("indexed " + count);
			}
		}
		this.table.finalizeObjectFile();
	}
	
	private Map indexSpecraByMass(List<Spectrum> specList){
		Map spectrumTable = new MultiValueMap();
		for(Iterator<Spectrum> it = specList.iterator(); it.hasNext();){
			Spectrum s = it.next();
			spectrumTable.put(getKey(s), s);
		}
		return spectrumTable;
	}
	

	
	private int getKey(Spectrum s){
		return (int)Math.round(s.parentMass);
	}
	
	private int getKey(double mass){
		return (int)Math.round(mass);
	}
	/**
	 * Retrieve a list of spectra by the parentmass of a query spectrum
	 * @param s
	 * @param tolerance
	 * @return
	 */
	public List<Spectrum> getSpectraByMass(Spectrum s, double tolerance){
		int key = getKey(s.parentMass);
		int width = (int)Math.ceil(tolerance);
		List<Spectrum> candidates = new ArrayList();
		for(int i = key-width; i <= key+width; i++){
			Object value = getSpectraByMass(i);
			if(value != null){
				candidates.addAll(getSpectraByMass(i));
			}
		}
//		this.avgTolerance = (19.0/20.0)*avgTolerance + (1/20.0)*tolerance; //keep track of a running tolerance for query
//		//System.out.println("avg tolerance: " + this.avgTolerance);
//		double ratio = Math.round(this.avgTolerance*2) / this.table.getKeepObjects();
//		if(ratio > 10 || ratio < 0.1){
//			this.table.setKeepObjects((int)Math.round(this.avgTolerance*2+1)*10);
//			//System.out.println("changing cache:  " + (int)Math.round(this.avgTolerance*2+1)*10);
//		}
		return candidates;
	}
	/**
	 * 
	 * @param mass
	 * @return
	 */
	public List<Spectrum> getSpectraByMass(double mass){
		int key = getKey(mass);
		List<Spectrum> cand = (List<Spectrum>)this.table.get(key);
		return cand;
	}
	
	/**
	 *This provide a convenient static way to index spectral library
	 * same as generateIndexedLibrary(String), using a CandidateSpectrumLibrary
	 * object
	 * @param libraryFile
	 * @param tolerance
	 * @return
	 */
	public static String indexSpectralLibrary(String libraryFile, double tolerance){
		String indexedFile = getIndexedLibraryName(libraryFile);
		if(!(new File(indexedFile).exists())){
			CandidateSpectrumGenerator gen = new CandidateSpectrumGenerator(libraryFile);
			gen.fragmentTolerance = tolerance;
			if(libraryFile.endsWith("mgf")){
				gen.generateIndexedLibrary(indexedFile);
			}else{
				gen.generateIndexedLibrary(indexedFile);
			}
		}else{
			System.out.println("indexed library file already existed, remove to reindexed");
		}
		return indexedFile;
	}
	
	/**
	 * Same as indexSpectralLibrary(String libraryFile, double tolerance), using 
	 * default tolerance (0.5 Da)
	 * @param libraryFile
	 * @return
	 */
	public static String indexSpectralLibrary(String libraryFile){
		return indexSpectralLibrary(libraryFile, 0.5);
	}
	
	/**
	 * 
	 * @param libraryFile
	 * @return
	 */
	public static String getIndexedLibraryName(String libraryFile){
		System.out.println("library file: " + libraryFile);
		String indexedFile = Utils.FileIOUtils.stripExtension(libraryFile) + ".map";
		System.out.println("generating indexed file: " + indexedFile);
		return indexedFile;
	}
	
	/**
	 * @param libraryFile
	 * @param indexedFile
	 */
	public static void indexSpectralLibrary(String libraryFile, String indexedFile){
		CandidateSpectrumGenerator gen = new CandidateSpectrumGenerator(libraryFile);
		gen.generateIndexedLibrary(indexedFile);
	}
	
	public static void testGenerator(){
		String libraryFile = "../mixture_linked/spectral_library/NIST_yeast_IT_v3.0_2009-05-04_7AA_decoy.sptxt";
		indexSpectralLibrary(libraryFile);
	}
	
	public static void testUseGenerator(){
		String libraryFile = "../mixture_linked/spectral_library/NIST_human_IT_2010-01-14_plusDecoy.mgf";
		CandidateSpectrumGenerator gen = new CandidateSpectrumGenerator(libraryFile);
		String libraryObjectFile = indexSpectralLibrary(libraryFile, 0.5);
		gen.loadLibraryFromFile(libraryObjectFile);
		List<Double> precursors = new ArrayList<Double>();
		for(int i = 0; i < 1000; i++){
			precursors.add(Math.random()*1500+500);
		}
		//Collections.sort(precursors);
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < precursors.size(); i++){
			Spectrum s = new Spectrum();
			s.parentMass = precursors.get(i);
			List<Spectrum> candidates = gen.getSpectraByMass(s, 3.0);
			System.out.println("number of candidates: "  + candidates.size());
		}
		System.out.println("matching " + 500 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void testStreamGenerator(){
		String libraryFile = "../mixture_linked/yeast0_annotated.mgf";
		indexSpectralLibrary(libraryFile);
	}
	

	public static void main(String[] args){
		//testGenerator();
		//testUseGenerator();
		//testStreamGenerator();
		//args[0] = "../mixture_linked/spectral_library/NIST_yeast_IT_v3.0_2009-05-04_7AA_decoy.sptxt";
		if(args.length != 2){
			System.out.println("java -Xmx3000M -cp MSPLIT.jar Spectrum.CandidateSpectrumGenerator <spectral library> <index File>");
			return;
		}
		String libraryFile = args[0];
		String indexFile = args[1];
		indexSpectralLibrary(libraryFile, indexFile);
	}
}	
