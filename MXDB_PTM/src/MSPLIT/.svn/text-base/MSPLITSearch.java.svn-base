package MSPLIT;

import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;

import MSPLIT.SpectrumIO.MGFParser;
import MSPLIT.SpectrumIO.MSPParser;

/**
 * Setup and run MSPLIT searches
 * @author bugouzhi
 *
 */
public class MSPLITSearch {
	private String libFile;
	private String[] queryFiles; //allow multiple files as query
	private double fragmentMassTolerance = 0.5;
	private Iterator<Spectrum> queriesIter;
	private List<Spectrum> library;
	private SpectrumLibSearcher searcher;
	public MSPLITSearch(String libFile, String[] queryFiles){
		this.libFile = libFile;
		this.queryFiles = queryFiles;
		System.out.println("Loading library: " + libFile);
		this.loadLibrary();
		System.out.println("Done loading library");
	}
	
	public MSPLITSearch(List<Spectrum> library, String[] queryFiles){
		this.libFile = "Default library-file name";
		this.queryFiles = queryFiles;
		this.library = library;
	}
	
	public MSPLITSearch(List<Spectrum> library, Iterator<Spectrum> iter){
		this.libFile = "Default library-file name";
		this.library = new ArrayList<Spectrum>();
		for(int i = 0; i < library.size(); i++){
			this.library.add(preprocessSpectrum(library.get(i)));
		}
		this.queriesIter = iter;
	}
	
	public double getFragmentMassTolerance() {
		return fragmentMassTolerance;
	}

	public void setFragmentMassTolerance(double fragmentMassTolerance) {
		this.fragmentMassTolerance = fragmentMassTolerance;
	}

	public List<Spectrum> getLibrary() {
		return library;
	}

	public void setLibrary(List<Spectrum> library) {
		this.library = library;
	}

	
	private void loadLibrary(){
		MSPParser parser = new MSPParser(this.libFile);
		this.library = new ArrayList<Spectrum>();
		while(parser.hasNext()){
			Spectrum s= parser.next();
			s = preprocessSpectrum(s);
			this.library.add(s);
		}
		System.out.println("loaded " + this.library.size() + " spectra from file");
	}
	
	private List<Spectrum> loadQueryFile(){
		MGFParser parser = new MGFParser(this.queryFiles[0]);
		List<Spectrum> queries = new ArrayList<Spectrum>();
		while(parser.hasNext()){
			Spectrum query = parser.next();
			queries.add(query);
		}
		this.queriesIter = queries.iterator();
		return queries;
	}
	
	private Spectrum preprocessSpectrum(Spectrum s){
		SpectrumUtil.sqrtSpectrum(s);
		SpectrumUtil.normalize(s);
		SpectrumUtil.scaleSpectrumMass(s, 0.9995);
		SpectrumUtil.toVector(s);
		return s;
	}
	
	public void search(){
		if(searcher == null){
			this.searcher = new SpectrumLibSearcher(this.library, 
				CosineSpectrumComparator.CosineComparator);
		}
		if(queriesIter == null){
			this.loadQueryFile();
		}
		long start = (new GregorianCalendar()).getTimeInMillis();
		while(this.queriesIter.hasNext()){
			Spectrum query = this.queriesIter.next();
			Spectrum query2 = query;//new SimpleSpectrum(query); //convert to vector
			preprocessSpectrum(query2);
			//searcher.topCandidate(query2);
			searcher.topCandidatePair(query2);
		}
		System.out.println("matching all spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMSPLITSearch(){
		String libFile = "..\\MSPLIB\\Lib\\human.msp";
		String queryFile = "..\\mixture_linked\\linked_peptide_spectra1.mgf";
		MSPLITSearch msplit = new MSPLITSearch(libFile, new String[]{queryFile});
		long start = (new GregorianCalendar()).getTimeInMillis();
		msplit.search();
		System.out.println("matching all spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testSimulatedMixture(){
		String libFile = "..\\MSPLIB\\Lib\\human.msp";
		MSPParser parser = new MSPParser(libFile);
		List specLib = new ArrayList();
		specLib.addAll(parser.readAllSpectra());
		System.out.println("read in library spectra: " + specLib.size());
		MixtureSpectrumGenerator generator = new MixtureSpectrumGenerator(specLib);
		System.out.println("start generating mixtures");
		List<Spectrum> queries = generator.generateRandomMixtures(1000, 0.1, 0.2);
		for(int i = 0; i < queries.size(); i++){
			AnnotatedMixtureSpectrum q = (AnnotatedMixtureSpectrum)queries.get(i);
			q.setSpectrumName(q.getPeptides()[0] + " & " + q.getPeptides()[1]);
		}
		System.out.println("generate query spectra: " + queries.size());
		MSPLITSearch msplit = new MSPLITSearch(specLib, queries.iterator());
		msplit.search();
	}
	
	public static void main(String[] args){
		//testMSPLITSearch();
		testSimulatedMixture();
	}
	
}
