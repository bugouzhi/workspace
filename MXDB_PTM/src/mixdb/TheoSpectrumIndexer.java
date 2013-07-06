package mixdb;
import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.Map;
import java.util.List;
import java.util.HashMap;
import java.util.Set;

import org.Spectrums.LargeHashMap;
import org.Spectrums.PeptideLite;
import org.Spectrums.Spectrum;

import SeqDB.DatabaseIndexer;
/**
 * implements the database interface, represent a collection of theoretical spectrum
 * @author Jian Wang
 *
 */
public class TheoSpectrumIndexer {
	private DatabaseIndexer pepDB;
	private LargeHashMap theoDB;
	private int minCharge = 2;
	private int maxCharge = 3;
	
	
	public TheoSpectrumIndexer(String dbPath){
		this.pepDB = new DatabaseIndexer(dbPath);
		createTheoSpectrum();
	}
	
	public TheoSpectrumIndexer(String dbPath, String objectPath){
		this.pepDB = new DatabaseIndexer(dbPath);
		this.theoDB = new LargeHashMap(objectPath);
		this.theoDB.loadLibraryFromFile(objectPath);
	}

	
	private void createTheoSpectrum(){
		this.theoDB = new LargeHashMap(pepDB.getDbPath()+".thmap");
		for(Iterator it = pepDB.getKeys().iterator(); it.hasNext();){
			Integer key = (Integer)it.next();
			List peps = pepDB.getPeptides(key, key, 1.0);
			List<Spectrum> theoList = new ArrayList<Spectrum>(peps.size());
			for(int i = 0; i < peps.size(); i++){
				PeptideLite peptide = (PeptideLite)peps.get(i);
				String p = pepDB.getSeq().getSubsequence(peptide.getBeginInd(), peptide.getEndInd()+1);
				for(int c = minCharge; c<=maxCharge; c++){
					Spectrum theo = TheoreticalSpectrumFactory.getTheoSpectrumX(p, c, TheoreticalSpectrumFactory.standardTypeMap, TheoreticalSpectrumFactory.standardIonMap);
					theoList.add(theo);
				}
			}
			theoDB.put(key, theoList);
		}
		this.theoDB.finalize();
	}
	
	public List<Spectrum> getTheoSpectrum(double fromMass, double toMass, double tolerance){
		int leftKey = pepDB.getKey(fromMass);
		int rightKey = pepDB.getKey(toMass);
		List<Spectrum> cands = new ArrayList<Spectrum>();
		for(int key = leftKey; key <= rightKey; key++){
			//System.out.println("key is: " + key);
			List<Spectrum> cand = (List<Spectrum>)this.theoDB.get(key);
			if(cand!=null){
				cands.addAll(cand);
			}
		}
		return cands;		
	}
	
	
	public static void testTheoDB(){
		String dbPath = "../mixture_linked/database/yeast_proteins_plusDecoy.fasta";
		String objectPath = "../mixture_linked/database/yeast_proteins_plusDecoy.fasta.thmap";
		TheoSpectrumIndexer theoIndexer = new TheoSpectrumIndexer(dbPath);
		//TheoSpectrumIndexer theoIndexer = new TheoSpectrumIndexer(dbPath, objectPath);
		Iterator it = theoIndexer.theoDB.getKeys().iterator();
		long start = (new GregorianCalendar()).getTimeInMillis();
		for(int i = 0; i < 1000 && it.hasNext(); i++){
			Integer key = (Integer)it.next();
			List values = theoIndexer.getTheoSpectrum(key, key, 1.0);
			System.out.println("At mass: " + key + " number of candidates: " + values.size());
			
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000);
	}
	
	public static void main(String[] args){
		testTheoDB();
	}
}
