package SeqDB;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jdbm.PrimaryHashMap;
import jdbm.PrimaryStoreMap;
import jdbm.RecordManager;
import jdbm.RecordManagerFactory;
import jdbm.SecondaryKeyExtractor;
import jdbm.SecondaryTreeMap;

import org.Spectrums.LargeHashMap;
import org.Spectrums.MZXMLReader;
import org.Spectrums.Mass;
import org.Spectrums.Spectrum;

import sequences.FastaSequence;

/**
 * Database indexer base on jdbm2
 * http://code.google.com/p/jdbm2/
 * rather than the self-wrote largehashmap :)
 * @author Jian Wang
 *
 */
public class Jdbm2Indexer {
	private static int MINPEPLENGTH=6;
	private static int MAXPEPLENGTH = 30;
	private String dbPath;
	private char[] ntermcut; //enzyme specificity, only allow simple rule, cut at specified residues
	private char[] ctermcut;
	private int numMissCut;
	private boolean ntermNonEnzy=false; //allow terminus to be non-enzymatic
	private boolean ctermNonEnzy=false;
	public double resolution=1.0; //this field should be set only during construction, ohterwise will mess up the key mapping	
	public static final char WILDCARD = '*';
	FastaSequence seq;
	PrimaryHashMap<Integer, List<PeptideLite>> peptideIndexes;
	PrimaryStoreMap<Long, PeptideLite> pepTable;
	SecondaryTreeMap<Integer, Long, PeptideLite> pepIndex;
	
	
	public Jdbm2Indexer(String dbPath){
		this(dbPath, 1.0);
	}
	
	public Jdbm2Indexer(String dbPath, double resolution){
		this(dbPath, new char[]{'R', 'K', '_'}, new char[]{'R', 'K'});
		//this(dbPath, new char[]{'F', 'Y', 'W', 'L', '_'}, new char[]{'F', 'Y', 'W', 'L'});
		this.resolution = resolution;
		try{
			indexDatabase();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public Jdbm2Indexer(String dbPath, char[] ntermcut, char[] ctermcut){
		this.dbPath = dbPath;
		this.ntermcut = ntermcut;
		this.ctermcut = ctermcut;
	}
	
	private boolean isIndexed(String dbPath){
		String indexObjectFile = getIndexFileName(dbPath);
		File f = new File(indexObjectFile);
		return f.isFile();
	}
	
	/**
	 * Generate a systematic name for the index file for the database
	 * we need to make sure this file name contain all the need information
	 * so the .map file can be used consistently
	 * information needed: 1) resolution
	 *                     2) enzyme spec
	 *                     3) length range
	 * @param dbPath
	 * @return
	 */
	private String getIndexFileName(String dbPath){
		StringBuffer indexFile = new StringBuffer(dbPath);
		indexFile.append("_len" + this.MINPEPLENGTH +"_" + this.MAXPEPLENGTH+"_spec_");
		for(int i = 0; i < this.ntermcut.length; i++){
			indexFile.append(this.ntermcut[i]);
		}
		indexFile.append("_");
		for(int i = 0; i < this.ntermcut.length; i++){
			indexFile.append(this.ntermcut[i]);
		}
		indexFile.append("_resolu_" + (int)(1000/this.resolution) + ".jbdm2map");
		return indexFile.toString();
	}
	
	public void indexDatabase() throws IOException{
		this.seq = new FastaSequence(this.dbPath);
		RecordManager recMan = RecordManagerFactory.createRecordManager(getIndexFileName(this.dbPath));
		this.pepTable = recMan.storeMap(getIndexFileName(this.dbPath));
		SecondaryTreeMap<Integer, Long, PeptideLite> pepInd = pepTable.secondaryTreeMap("pmassIndex",
				new SecondaryKeyExtractor<Integer, Long, PeptideLite>(){
						public Integer extractSecondaryKey(Long key, PeptideLite p){
							return p.getParentMass();
						}
				});

		if(isIndexed(this.dbPath)){
			return; // no need to create index
		}
		//this.peptideIndexes = table;
		long beginIndex = 1;
		double currentMass = 0.0;
		int count = 0;
		for(long size = this.seq.getSize(); beginIndex < size-MAXPEPLENGTH;){
			//System.out.println(beginIndex);
			currentMass = 0.0;
			//System.out.println("first cut: " + seq.getCharAt(beginIndex-1));
			for(long j = 0; j < MAXPEPLENGTH; j++){
				char c = seq.getCharAt(beginIndex+j);
				//System.out.println("current char: " + c);
				currentMass += Mass.getAAMass(c);
				if(seq.isTerminator(beginIndex+j)){
					currentMass -= Mass.getAAMass(c);
					int key = getKey(currentMass);
					pepTable.putValue(new PeptideLite((int)beginIndex, (int)(beginIndex+j-1), key));
					count++;
					break;
				}
				
				if(j < MINPEPLENGTH-1 || currentMass > 10000){
					continue;
				}
				
				if(checkCterm(c)){
					//System.out.println("second cut: " + seq.getCharAt(beginIndex+j));
					int key = getKey(currentMass);
					//System.out.println(this.seq.getSubsequence(beginIndex, beginIndex+j+1) + "\t" + currentMass);
					//System.out.println("current mass is: " + currentMass);
					pepTable.putValue(new PeptideLite((int)beginIndex, (int)(beginIndex+j), key));
					count++;
				}
				if(count % 100000 == 0){
					recMan.commit();
					System.out.println("Processed: " + count);
					//recMan.close();
					//return;
				}
			}
			beginIndex = nextBeginIndex(beginIndex);
		}
		System.out.println("Done indexing, indexed peptides: " + count);
		recMan.commit();
		//recMan.close();
	}
	
	public Collection getKeys(){
		return this.peptideIndexes.keySet();
	}
	
	public List<PeptideLite> getPeptides(double fromMass, double toMass, double tolerance){
		//System.out.println("from: " + fromMass + " to " + toMass);
		return getPeptidesWithC13(fromMass, toMass, tolerance, 0);
			
	}
	
	public List<PeptideLite> getPeptidesWithC13(double fromMass, double toMass, double tolerance, int C13){
		//System.out.println("from: " + fromMass + " to " + toMass);
		List<PeptideLite> cands = new ArrayList<PeptideLite>();
		double offset = Mass.C13 - Mass.C12;
		if(tolerance > offset){   //no need to consider this when tolerance is large
			C13 = 0;
		}
		for(int i = 0; i <= C13; i++){
			double leftMass = fromMass - offset*C13;
			double rightMass = toMass - offset*C13;
			int leftKey = getKey(leftMass);
			int rightKey = getKey(rightMass);
			for(int key = leftKey; key <= rightKey; key++){
				//System.out.println("key i- s: " + key);
				List<PeptideLite> cand = (List<PeptideLite>)this.peptideIndexes.get(key);
				for(PeptideLite pep:this.pepIndex.getPrimaryValues(key)){
					cands.add(pep);
				}
			}
		}
		//System.out.println("found: " + cands.size());
		return cands;
			
	}
	
	public FastaSequence getSeq() {
		return seq;
	}

	public void setSeq(FastaSequence seq) {
		this.seq = seq;
	}

	private long nextBeginIndex(long start){
		for(long i = start, size = this.seq.getSize(); i < size; i++){
			for(int j = 0; j < ntermcut.length; j++){
				if(this.seq.getCharAt(i) == ntermcut[j]){
					return i+1;
				}
			}
		}
		return this.seq.getSize();
	}
	
	private boolean checkCterm(char c){
		for(int j = 0; j < ctermcut.length; j++){
			if(c == ctermcut[j]){
				return true;
			}
		}
		return false;
	}
	
	public int getKey(double mass){
		return (int)Math.round(mass/resolution)+1; 
	}
	
	
	
	public String getDbPath() {
		return dbPath;
	}

	public void setDbPath(String dbPath) {
		this.dbPath = dbPath;
	}
	
	
	public static void testIndexDB(){
		double bin = 0.01;
		Jdbm2Indexer indexer = new Jdbm2Indexer("../mixture_linked/database/Human_allproteins_plusDecoy.fasta", bin);
		String specFile = "..//mixture_linked/yeast_data/klc_010908p_yeast-digest.mzXML";
		MZXMLReader reader = new MZXMLReader(specFile);
		//for(int i = 0; i < 1000; i++){
		for(;reader.hasNext();){
			Spectrum s = reader.next();
			//double precursor = 600 + 2000*Math.random();
			//int charge = (int)(2 + 3*Math.random());
			//precursor = precursor /charge;
			int i = s.scanNumber;
			double precursor = s.parentMass;
			int charge = s.charge;
			double mass =  PeptideUtils.getPeptideMass(precursor, charge);
			double ppm = 30;
			double tolerance = mass*ppm/1000000;
			//System.out.println("tolerance: " + tolerance);
			List<PeptideLite> cands = indexer.getPeptidesWithC13(mass - tolerance, 
					mass + tolerance, tolerance, 1);
			System.out.println("Spectrum: " + i + "\t" + precursor + "\t" + charge +"\tpeptides: " + cands.size());
		}
	}
	
	
	public static void main(String[] args){
		testIndexDB();
	}
	
}
