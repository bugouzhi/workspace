package org.Spectrums;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Create decovuluted spectrum
 * @author Jian Wang
 *
 */
public class DeconvolutedSpectrum extends Spectrum{
	Spectrum orig; //the orignal spectrum
	double fragmentTolerance = 30;
	
	public DeconvolutedSpectrum(Spectrum s){
		super(s);
		this.orig = s;
		deconvolute();
	}
	
	private void deconvolute(){
		int[] peakCharge = this.compuatePeakCharge();
		List<Peak> pList = this.getPeak();
		double isomass = Mass.C13 - Mass.C12;
		for(int i = 0; i < pList.size()-1; i++){
			Peak current = pList.get(i);
			System.out.println("current: " + current + " peak charge: " + peakCharge[i]);
			
			if(peakCharge[i] > 1){
				current.setMoz(current.getMass()*peakCharge[i] - (peakCharge[i]-1)*isomass);
			}
		}
		Collections.sort(pList, PeakMassComparator.comparator);
		this.setPeaks(pList);
		removeIsoPeaks();
	}
	
	private int[] compuatePeakCharge(){
		int[] peakCharge = new int[this.getPeak().size()];
		List<Peak> pList = this.getPeak();
		for(int i = 0; i < pList.size()-1; i++){
			Peak current = pList.get(i);
			int j = i+1;
			Peak next = pList.get(j);
			double diff = next.getMass() - current.getMass();
			while(diff < 1.01){
				for(int charge=this.orig.charge; charge > 0; charge--){
					if(Math.abs(diff - Mass.PROTON_MASS/charge)*1000000/ current.getMass() < this.fragmentTolerance){
						peakCharge[i] = charge;
						break;
					}
				}
				j++;
				next = pList.get(j);
				diff = next.getMass() - current.getMass();
			}
		}
		return peakCharge;
	}
	
	private void removeIsoPeaks(){
		List<Peak> pList = this.getPeak();
		List<Peak> toBeRemoved = new ArrayList<Peak>();
		for(int i = 0; i < pList.size()-1; i++){
			Peak current = pList.get(i);
			//System.out.println("current peak: " + current);
			int j = i+1;
			Peak next = pList.get(j);
			double diff = next.getMass() - current.getMass();
			int isocount = 1;
			while(diff < 1.01 && j < pList.size()-2){
				for(int charge=this.orig.charge; charge > 0; charge--){
					if(Math.abs(diff - isocount*Mass.PROTON_MASS/charge)*1000000/ current.getMass() < this.fragmentTolerance){
						current.setIntensity(current.getIntensity()+next.getIntensity());
						//System.out.println("removing iso peaks @ " + charge + ":\t" +  next);		
						break;
					}
				}
				toBeRemoved.add(next);
				isocount++;
				j++;
				next = pList.get(j);
				diff = next.getMass() - current.getMass();
				
			}
		}
		//System.out.println("before: " + pList.size());
		//System.out.println("removed peaks: " + toBeRemoved);
		pList.removeAll(toBeRemoved);
		//System.out.println("after: " + pList.size());
	}
	
	public static void testDecolvolution(){
		String spectrumFile = "..\\mixture_linked/linked_peptide_library/sumo_lib/Veronica_HCD/20111122_ananiav_Sumo_training_library_0,5ul_lib2_HCD40.mzXML";
		int scan = 3910;
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		Spectrum s = reader.getSpectrum(scan);
		DeconvolutedSpectrum Ds = new DeconvolutedSpectrum(s);	
	}
	
	public static void main(String[] args){
		testDecolvolution();
	}
	
}
