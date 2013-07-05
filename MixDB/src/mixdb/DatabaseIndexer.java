package mixdb;
import org.Spectrums.LargeHashMap;
import org.Spectrums.PeptideLite;
import org.Spectrums.Mass;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.List;
import java.util.HashMap;
import sequences.*;
/**
 * A simple database indexer for protein sequences, index all peptides
 * @author Jian Wang
 *
 */
public class DatabaseIndexer {
	private static int MINPEPLENGTH=6;
	private static int MAXPEPLENGTH = 30;
	private String dbPath;
	private char[] ntermcut; //enzyme specificity, only allow simple rule, cut at specified residues
	private char[] ctermcut;
	private int numMissCut;
	private boolean ntermNonEnzy=false; //allow terminus to be non-enzymatic
	private boolean ctermNonEnzy=false;
	private double resolution=1.0; //this field should be set only during construction, otherwise will mess up the key mapping	
	FastaSequence seq;
	LargeHashMap peptidesIndex;
	
	public DatabaseIndexer(String dbPath){
		this(dbPath, new char[]{'R', 'K', '_'}, new char[]{'R', 'K'});
		indexDatabase();
	}
	
	public DatabaseIndexer(String dbPath, char[] ntermcut, char[] ctermcut){
		this.dbPath = dbPath;
		this.ntermcut = ntermcut;
		this.ctermcut = ctermcut;
	}
	
	public void indexDatabase(){
		this.seq = new FastaSequence(this.dbPath);
		Map<Integer, List<PeptideLite>> table = new HashMap();
		long beginIndex = 1;
		double currentMass = 0.0;
		int count = 0;
		File index = new File(this.dbPath+".map");
		if(index.exists()){
			this.peptidesIndex = new LargeHashMap(this.dbPath+".map");
			this.peptidesIndex.loadLibraryFromFile(this.peptidesIndex.getLibraryObjectFile());
			return;
		}
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
					List<PeptideLite> candList;
					if(table.containsKey(key)){
						candList=table.get(key);
					}else{
						candList = new ArrayList();
					}
					//System.out.println("end: " + this.seq.getSubsequence(beginIndex, beginIndex+j) + "\t" + currentMass);
					candList.add(new PeptideLite((int)beginIndex, (int)(beginIndex+j-1)));
					table.put(key, candList);
					count++;
					break;
				}
				
				if(j < MINPEPLENGTH-1 || currentMass > 10000){
					continue;
				}
				
				if(checkCterm(c)){
					//System.out.println("second cut: " + seq.getCharAt(beginIndex+j));
					int key = getKey(currentMass);
					List<PeptideLite> candList;
					if(table.containsKey(key)){
						candList=table.get(key);
					}else{
						candList = new ArrayList();
					}
					//System.out.println(this.seq.getSubsequence(beginIndex, beginIndex+j+1) + "\t" + currentMass);
					//System.out.println("current mass is: " + currentMass);
					candList.add(new PeptideLite((int)beginIndex, (int)(beginIndex+j)));
					table.put(key, candList);
					count++;
				}
			}
			beginIndex = nextBeginIndex(beginIndex);
		}
		System.out.println("Done indexing, indexed peptides: " + count);
		this.peptidesIndex = new LargeHashMap(this.dbPath+".map");
		this.peptidesIndex.buildTable(table);
		this.peptidesIndex.loadLibraryFromFile(this.peptidesIndex.getLibraryObjectFile());
	}
	
	public Collection getKeys(){
		return this.peptidesIndex.getKeys();
	}
	
	public List<PeptideLite> getPeptides(double fromMass, double toMass, double tolerance){
		int leftKey = getKey(fromMass);
		int rightKey = getKey(toMass);
		List<PeptideLite> cands = new ArrayList<PeptideLite>();
		for(int key = leftKey; key <= rightKey; key++){
			//System.out.println("key is: " + key);
			List<PeptideLite> cand = (List<PeptideLite>)this.peptidesIndex.get(key);
			if(cand!=null){
				cands.addAll(cand);
			}
		}
		for(int i = 0; i < cands.size(); i++){
			PeptideLite pep = cands.get(i);
			pep.setFastaseq(this.getSeq());
		}
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
		DatabaseIndexer indexer = new DatabaseIndexer("../mixture_linked/database/Human_allproteins_plusDecoy.fasta");
	}
	
	
	public static void main(String[] args){
		testIndexDB();
	}
	
	
	
	

}
