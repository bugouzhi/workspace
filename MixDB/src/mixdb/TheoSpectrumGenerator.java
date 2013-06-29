package mixdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;
import org.Spectrums.Mass;
import org.Spectrums.Peptide;
import org.Spectrums.PeptideLite;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumMassComparator;

import sequences.FastaSequence;

/**
 * A class that takes in a peptide database, it generate theoretical spectrum
 * for the peptides candidates, to be efficient it store a subsset of the theoretical
 * spectrum in memory, so it won't need to re-create the theo-spectrum repeatedly
 * @author Jian Wang
 *
 */
public class TheoSpectrumGenerator {
	DatabaseIndexer pepDB;
	TreeSet<Spectrum> candidates;
	private double currentMinMZ;
	private double currentMaxMZ;
	private int minCharge = 2;
	private int maxCharge = 3;
	private boolean matchCharge = false;
	private boolean DEBUG = false;
	public boolean isMatchCharge() {
		return matchCharge;
	}

	public void setMatchCharge(boolean matchCharge) {
		this.matchCharge = matchCharge;
	}

	public TheoSpectrumGenerator(DatabaseIndexer pepDB){
		this.pepDB = pepDB;
		init();
	}
	
	
	private void init(){
		this.candidates = new TreeSet<Spectrum>(SpectrumMassComparator.comparator);
	}
	
	
	public Collection<Spectrum> getCandidates(Spectrum s, double parentMassTolerance){
		List<Spectrum> theoCands = new ArrayList();
		for(int c = minCharge; c <= maxCharge; c++){
			double pm = s.parentMass*c-Mass.WATER-c*Mass.PROTON_MASS;
			double tolerance = parentMassTolerance*c;
			double currentLeft = currentMaxMZ*c-Mass.PROTON_MASS*c-Mass.WATER;
			currentLeft = pm-tolerance > currentLeft ? pm-tolerance : currentLeft;
			double currentRight = pm+tolerance;
			List<PeptideLite> pep = null;
			if(DEBUG){
				System.out.println("left is: " + currentLeft + " right is: " + currentRight);
			}
			if(currentLeft < currentRight){
				pep = this.pepDB.getPeptides(currentLeft, currentRight, tolerance);
			}
			if(pep != null){
				theoCands.addAll(TheoreticalSpectrumFactory.getArrayTheoreticalSpectrum(pep, this.pepDB.getSeq(), s, c));
			}
		}
		//remove candidates
		Spectrum current=null;
		if(DEBUG){
			System.out.println("number of candidate before: " + this.candidates.size());
		}
		for(Iterator<Spectrum> iter = this.candidates.iterator(); iter.hasNext();){
			current = (Spectrum)iter.next();
			if(Math.abs(s.parentMass - current.parentMass) > parentMassTolerance 
					|| (matchCharge && s.charge != current.charge)){
				iter.remove();
			}
		}
		if(DEBUG){
			System.out.println("number of candidate after: " + this.candidates.size());
			System.out.println("current set: [ " + this.currentMinMZ + " , " + this.currentMaxMZ +" ] ");
		}
		//this.candidates.clear();
		int added=0;
		for(Iterator<Spectrum> iter = theoCands.iterator(); iter.hasNext();){
			Spectrum cand = iter.next();
			if(Math.abs(cand.parentMass - s.parentMass) <= parentMassTolerance
					&& cand.parentMass > this.currentMaxMZ
					&& (!matchCharge || s.charge == current.charge)){
				this.candidates.add(cand);
				added++;
			}	
		}
		if(DEBUG){
			System.out.println("added candidates: " + added);
			System.out.println("candidates has size: " + this.candidates.size());
		}
		if(this.candidates.size() > 0){
			this.currentMinMZ = this.candidates.first().parentMass;
			this.currentMaxMZ = this.candidates.last().parentMass;
			if(DEBUG){
				System.out.println("current set: [ " + this.currentMinMZ + " , " + this.currentMaxMZ +" ] ");
			}
		}
		return this.candidates;
	}
	
	/**
	 * This version of the method handle querying multiple precursors with same spectrum
	 * the windowWidth tells how wide the isolation window for selecting precursors for fragmentation
	 * for ms/ms 
	 * @param s
	 * @param parentMassTolerance
	 * @param windowWidth
	 * @return
	 */
	public Collection<Spectrum> getCandidates(Spectrum s, double windowWidth, double precursorMass, double tolerance, int charge){
		Collection<Spectrum> wideCandidates = this.getCandidates(s, windowWidth);
		System.out.println("window candidates: " + wideCandidates.size());
		List<Spectrum> filteredCandidates = new ArrayList<Spectrum>();
		for(Iterator<Spectrum> iter = wideCandidates.iterator(); iter.hasNext();){
			Spectrum cand = iter.next();
			if(Math.abs((cand.parentMass - precursorMass)) < tolerance
					&& cand.charge == charge){
				filteredCandidates.add(cand);
			}
		}
		return filteredCandidates;
	}
	
	
	
	public static void testTheoGenerator(){
		
	}
	
	public static void main(String[] args){
		
	}
	
}
