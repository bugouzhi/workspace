package org.Spectrums;

import java.util.List;
import java.util.Set;
import java.util.Iterator;
import java.util.Map;
import java.util.HashMap;
import org.jgrapht.graph.SimpleGraph;

import Utils.ArrayUtils;

import java.util.ArrayList;
import java.util.Collections;
/**
 * Analyze annotated spectrum, extract ions type information or statistics
 * @author jian wang
 *
 */
public class SpectrumIonsLearner {
	public static int MAX_CHARGE_VALUE = 6;
	private SpectrumLib lib;
	public SpectrumIonsLearner(SpectrumLib lib){
		this.lib = lib;
	}
	
	public void getIonsStatistics(){
		//lib.windowFilterPeaks(6, 25);
		int maxCharge = 5;
		for(int charge = 1; charge <= maxCharge; charge++){
			System.out.println("Analyzing peptides with charge: " +   charge);
			getIonsStatistics2(charge);
		}
	}
	
	public void getNoisesStatistics(){
		//lib.windowFilterPeaks(6, 25);
		int maxCharge = 5;
		for(int charge = 1; charge <= maxCharge; charge++){
			System.out.println("Analyzing peptides with charge: " +   charge);
			getNoiseStatistics(charge);
		}
	}

	public int getIonsStatistics(int pepCharge){
		List list = lib.getAllSpectrums();
		Spectrum s;
		Peptide p;
		double[] fraction;
		double[] fractions = new double[Mass.standardIonsType.length*pepCharge];
		Map counts = null;
		int count = 0;
		for(int i = 0; i < list.size(); i++){
			s = (Spectrum)list.get(i);
			if(s.peptide.contains("Scan") || s.charge != pepCharge){
				continue;
			}
			fraction = getIonsFraction(s, pepCharge);
			ArrayUtils.addArray(fractions, fraction);
			count++;
			//return;
		}
		ArrayUtils.divideArray(fractions, (double)count);
		System.out.println("Total spectra gathered: " + count);
		System.out.println("Ion-Type: " + "averge frequencey");
		for(int c = 0; c < pepCharge; c++){
			for(int i = 0 ; i < Mass.standardIonsType.length; i++){
				System.out.println(Mass.standardIonsType[i] + "@" + (c+1)+ "@"+pepCharge + "\t" + fractions[c*Mass.standardIonsType.length+i]);
			}
		}
		System.out.println();
		return count;
	}
	
	public int getIonsStatistics2(int pepCharge){
		List list = lib.getAllSpectrums();
		Spectrum s;
		Peptide p;
		double[] fraction;
		double[] fractions = new double[Mass.standardIonsType.length*pepCharge];
		Map counts = new HashMap();
		int total = 0;
		int count = 0;
		for(int i = 0; i < list.size(); i++){
			s = (Spectrum)list.get(i);
			if(s.peptide.contains("Scan") || s.charge != pepCharge){
				continue;
			}
			TheoreticalSpectrum t = new TheoreticalSpectrum(new Peptide(s.peptide));
			SimpleGraph matchingG = t.matchSpectrum2(s, 0.5);
			this.getIonsCount(matchingG, counts);
			count++;
			total += s.peptide.split("\\.")[0].length();
			//return;
		}
		for(int c = 1; c <= pepCharge; c++){
			for(int i = 0; i<Mass.standardIonsType.length; i++){
				if(counts.containsKey(Mass.standardIonsType[i]+"@"+c)){
					fractions[(c-1)*Mass.standardIonsType.length+i] = ((Integer)counts.get(Mass.standardIonsType[i]+"@"+c)).intValue();
					//System.out.println("count is: " + fractions[i] + " peptide is: " + p.getPeptide());
				}
			}
		}
		ArrayUtils.divideArray(fractions, (double)total);
		System.out.println("Total spectra gathered: " + count);
		System.out.println("Ion-Type: " + "averge frequencey");
		for(int c = 0; c < pepCharge; c++){
			for(int i = 0 ; i < Mass.standardIonsType.length; i++){
				System.out.println(Mass.standardIonsType[i] + "@" + (c+1)+ "@"+pepCharge + "\t" + fractions[c*Mass.standardIonsType.length+i]);
			}
		}
		System.out.println();
		return count;
	}
	
	public void getNoiseStatistics(int pepCharge){
		List list = lib.getAllSpectrums();
		Spectrum s;
		TheoreticalSpectrum t;
		SimpleGraph matchingG;
		int count = 0;
		Iterator it;
		Peak p;
		Map<Integer, Integer>rankCounts = new HashMap<Integer, Integer>();
		for(int i = 0; i < list.size(); i++){
			List<Peak> noisePeaks = new ArrayList<Peak>();
			s = (Spectrum)list.get(i);
			s.computePeakRank();
			if(s.peptide.contains("Scan") || s.charge != pepCharge){
				continue;
			}
			t = new TheoreticalSpectrum(new Peptide(s.peptide));
			matchingG = t.matchSpectrum2(s, 0.5);
			it = matchingG.vertexSet().iterator();
			while(it.hasNext()){
				p = (Peak)it.next();
				if(!(p instanceof LabelledPeak) && matchingG.degreeOf(p) == 0){
					noisePeaks.add(p);
				}			
			}
			
			for(int j = 0, size2 = noisePeaks.size(); j< size2; j++){
				//System.out.println(i+s.peptide + " noise peak has rank " + noisePeaks.get(j) + "\t" + noisePeaks.get(j).getRank() + " total: " + s.getPeak().size());
				Integer rank = new Integer(SpectrumUtil.normalizedRank(noisePeaks.get(j), s));
				if(rankCounts.containsKey(rank)){
					Integer rankcount = rankCounts.get(rank);
					rankcount++;
					rankCounts.put(rank, rankcount);
				}else{
					rankCounts.put(rank, new Integer(1));
				}
			}
			count++;
		}
		List<Integer> keys = new ArrayList(rankCounts.keySet().size());
		keys.addAll(rankCounts.keySet());
		Collections.sort(keys);
		for(Iterator<Integer> iter = keys.iterator(); iter.hasNext();){
			Integer rank = iter.next();
			System.out.println("peak rank:\t" + rank + "\thas freq:\t" 
					+ ((rankCounts.get(rank)).doubleValue()/(double)count));
		}
		
		//return;
	}
	
	/**
	 * Noise model with intensity base rather than rank base key
	 */
	public void getNoiseStatistics2(){
		List list = lib.getAllSpectrums();
		Spectrum s;
		TheoreticalSpectrum t;
		SimpleGraph matchingG;
		int count = 0;
		Iterator it;
		Peak p;
		Map<Integer, Integer>rankCounts = new HashMap<Integer, Integer>();
		for(int i = 0; i < list.size(); i++){
			List<Peak> noisePeaks = new ArrayList<Peak>();
			s = (Spectrum)list.get(i);
			s.computePeakRank();
			t = new TheoreticalSpectrum(new Peptide(s.peptide));
			matchingG = t.matchSpectrum2(s, 0.5);
			it = matchingG.vertexSet().iterator();
			while(it.hasNext()){
				p = (Peak)it.next();
				if(!(p instanceof LabelledPeak) && matchingG.degreeOf(p) == 0){
					noisePeaks.add(p);
				}			
			}
			
			for(int j = 0, size2 = noisePeaks.size(); j< size2; j++){
				//System.out.println(i+s.peptide + " noise peak has rank " + noisePeaks.get(j) + "\t" + noisePeaks.get(j).getRank() + " total: " + s.getPeak().size());
				Integer rank = new Integer(SpectrumUtil.normalizedRank(noisePeaks.get(j), s));
				if(rankCounts.containsKey(rank)){
					Integer rankcount = rankCounts.get(rank);
					rankcount++;
					rankCounts.put(rank, rankcount);
				}else{
					rankCounts.put(rank, new Integer(1));
				}
			}
			count++;
		}
		List<Integer> keys = new ArrayList(rankCounts.keySet().size());
		keys.addAll(rankCounts.keySet());
		Collections.sort(keys);
		for(Iterator<Integer> iter = keys.iterator(); iter.hasNext();){
			Integer rank = iter.next();
			System.out.println("peak rank:\t" + rank + "\thas freq:\t" 
					+ ((rankCounts.get(rank)).doubleValue()/(double)count));
		}
		
		//return;
	}
	
	public double[] getIonsFraction(Spectrum s, int pepCharge){
		Map counts;
		TheoreticalSpectrum th;
		SimpleGraph matchingG;
		Peptide p = new Peptide(s.peptide);
		th = new TheoreticalSpectrum(p);
		matchingG = th.matchSpectrum2(s, 0.5);
		counts = getIonsCount(matchingG);
		double[] fractions = new double[Mass.standardIonsType.length*pepCharge];
		for(int c = 1; c <= pepCharge; c++){
			for(int i = 0; i<Mass.standardIonsType.length; i++){
				if(counts.containsKey(Mass.standardIonsType[i]+"@"+c)){
					fractions[(c-1)*Mass.standardIonsType.length+i] = ((Integer)counts.get(Mass.standardIonsType[i]+"@"+c)).intValue();
					//System.out.println("count is: " + fractions[i] + " peptide is: " + p.getPeptide());
					fractions[(c-1)*Mass.standardIonsType.length+i] /= (double) (p.getPeptide().length());
					//System.out.println("count is: " + fractions[i] + " peptide is: " + p.getPeptide());
				}
			}
		}
		return fractions;
	}
	
	private Map getIonsCount(SimpleGraph g){
		HashMap<String, Integer> counts = new HashMap();
		Set vertices = g.vertexSet();
		Iterator it = vertices.iterator();
		Peak p;
		LabelledPeak lp;
		Integer c, one = new Integer(1);
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (LabelledPeak)p;
				if(g.degreeOf(p) > 0){
					if(counts.containsKey(lp.getType()+"@"+lp.getCharge())){
						c=(Integer)counts.get(lp.getType()+"@"+lp.getCharge());
						c++;
						counts.put(lp.getType()+"@"+lp.getCharge(), c);
					}else{
						counts.put(lp.getType()+"@"+lp.getCharge(), one);
					}
				}
			}			
		}
		return counts;
	}
	
	private Map getIonsCount(SimpleGraph g, Map counts){
		Set vertices = g.vertexSet();
		Iterator it = vertices.iterator();
		Peak p;
		LabelledPeak lp;
		Integer c, one = new Integer(1);
		while(it.hasNext()){
			p = (Peak)it.next();
			if(p instanceof LabelledPeak){
				lp = (LabelledPeak)p;
				if(g.degreeOf(p) > 0){
					if(counts.containsKey(lp.getType()+"@"+lp.getCharge())){
						c=(Integer)counts.get(lp.getType()+"@"+lp.getCharge());
						c++;
						counts.put(lp.getType()+"@"+lp.getCharge(), c);
					}else{
						counts.put(lp.getType()+"@"+lp.getCharge(), one);
					}
				}
			}			
		}
		return counts;
	}
		
	public static void testGetIonStat(){
//		String file = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
		String file = ".\\MSPLib\\Lib\\ecoli.msp";
		String annotation = ".\\mixture_linked\\trps\\result.txt";
//		SpectrumLib lib = new SpectrumLib(file, "MGF");
		SpectrumLib lib = new SpectrumLib(file, "MSP");
		lib.removeModSpectra();
		lib.annoateSpectrumFromInspectFile(annotation);
		SpectrumIonsLearner learner = new SpectrumIonsLearner(lib);
		learner.getIonsStatistics();
		//learner.getNoisesStatistics();
	}
	
	public static void main(String[] args){
		testGetIonStat();
	}
}
