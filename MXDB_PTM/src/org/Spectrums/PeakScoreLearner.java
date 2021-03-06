package org.Spectrums;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.graph.SimpleGraph;

/**
 * Similar to SpectrumIonsLearner, learn how to score a matching peaks
 * back by a generic lookup table
 * @author Jian Wang
 *
 */
public abstract class PeakScoreLearner {
	private String[] ionsType = Mass.standardIonsType;
	private HashMap<String, Integer> ionIndex;
	private int minCharge = 1;
	private int maxCharge = 4;
	private LookUpTable table;
	private SpectrumLib annotatedSet;
	public PeakScoreLearner(SpectrumLib lib){
		this.annotatedSet = lib;
		this.initializeIonIndexTable();
		this.table = initializeTable();
	}
	
	public PeakScoreLearner(SpectrumLib lib, String[] ionsType, int minCharge, int maxCharge){
		this.annotatedSet = lib;
		this.ionsType = ionsType;
		this.minCharge = minCharge;
		this.maxCharge = maxCharge;
		this.initializeIonIndexTable();
		this.table = initializeTable();
	}
	
	/**
	 * One of the two method  a learner needs to implement,
	 * in order to specify what kind of model is desired
	 */
	private LookUpTable initializeTable(){
		LookUpTable table = new LookUpTable( 
			new int[] {maxCharge,  maxCharge, this.ionsType.length});
		return table;
	}
	
	private void initializeIonIndexTable(){
		this.ionIndex = new HashMap<String, Integer>();
		for(int i = 0; i < ionsType.length; i++){
			ionIndex.put(ionsType[i], new Integer(i));
			//System.out.println("storing ion: " + ionsType[i]);
		}
	}
	
	public int getIonIndex(LabelledPeak lp){
		if(!this.ionIndex.containsKey(lp.getType())){
			throw new IllegalArgumentException("Invalide ion type " + lp.getType());
		}
		return ionIndex.get(lp.getType()).intValue();
	}
	
	/**
	 * One of the two method a learner needs to implement,
	 * in order to specify what kind of model is desired
	 * @param lp
	 * @return
	 */
	
	public int[] getIndex(LabelledPeak lp){
		int peptideCharge = lp.getPep().getCharge()-1;
		int peakCharge = lp.getCharge()-1;
		int ionIndex = getIonIndex(lp);
		return new int[]{peptideCharge, peakCharge, ionIndex};
	}
	
	public double getValue(int[] index){
		return this.table.get(index);
	}
	
	public void getIonsCount(){
		List list = this.annotatedSet.getAllSpectrums();
		Spectrum s;
		Peptide p;
		int total = 0;
		int count = 0;
		LookUpTable totalCount = initializeTable();
		for(int i = 0; i < list.size(); i++){
			s = (Spectrum)list.get(i);
			//System.out.println("peptide is: " + s.peptide);
			TheoreticalSpectrum t = new TheoreticalSpectrum(new Peptide(s.peptide));
			SimpleGraph matchingG = t.matchSpectrum2(s, 0.5);
			this.getIonsCount(matchingG, totalCount);
			//return;
		}
		LookUpTable.mapOperator(this.table, totalCount, TableElementOperator.Divider.d);
		//LookUpTable.mapOperator(this.table, 0.05, TableElementOperator.Divider.d);
		//LookUpTable.mapOperator(this.table, 1, TableElementOperator.Log.l);
		printIonTable();
	}
	
	
	private void getIonsCount(SimpleGraph g, LookUpTable totalCount){
		Set vertices = g.vertexSet();
		Iterator it = vertices.iterator();
		Peak p;
		LabelledPeak lp;
		Integer c, one = new Integer(1);
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (LabelledPeak)p;
				int[] index = getIndex(lp);
				if(g.degreeOf(p) > 0){
					this.table.incrementIonCount(index);
				}
				totalCount.incrementIonCount(index);
			}
		}
	}
	
	private void printIonTable(){
		for(int pepCharge = 0; pepCharge < 3; pepCharge++){
			for(int peakCharge = 0; peakCharge < 3; peakCharge++){
				for(int ionIndex = 0; ionIndex < this.ionsType.length; ionIndex++){
					int[] index = {pepCharge, peakCharge, ionIndex};
					System.out.println(this.ionsType[ionIndex] + "@" + (peakCharge+1) + "@" + (pepCharge+1) + ": " + this.table.get(index));
				}
			}
		}
	}
	
	public static void testGetIonStat(){
		String file = ".\\MSPLib\\Lib\\ecoli.msp";
		SpectrumLib lib = new SpectrumLib(file, "MSP");
		lib.removeModSpectra();
		//PeakScoreLearner learner = new PeakScoreLearner(lib);
		//learner.getIonsCount();
		//learner.getNoisesStatistics();
	}
	
	public static void main(String[] args){
		testGetIonStat();
	}
	
}
