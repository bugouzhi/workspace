package org.Spectrums;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

import mixdb.ArrayTheoreticalSpectrum;
import mixdb.LinkedTheoSpectrumFactory;

//lazy-evaluated version of theoretical spectrum
//only generate the spectrum upon requested
//this is so to handle large spectrum library
public class LazyEvaluateLinkedSpectrum extends TheoreticalSpectrum{
	private int linkedCharge;
	Peptide p1;
	Peptide p2;
	double linkerMass;
	public int getLinkedCharge() {
		return linkedCharge;
	}
	
	public void setLinkedCharge(int linkedCharge) {
		this.linkedCharge = linkedCharge;
	}
	
	
	//think about what is proper way to do this? is it to initialize everything but not created the peasks?
	public LazyEvaluateLinkedSpectrum(Peptide p, int linkedCharge){
		this.linkedCharge = linkedCharge;
		this.peptide = p.toString();//p.getPeptide()+"."+p.getCharge();
		this.p = p;
		this.p1 = p;
		this.charge = linkedCharge;
		this.parentMass = ((p.getParentmass()-p.getCharge()*Mass.PROTON_MASS)+linkedCharge*Mass.PROTON_MASS)/linkedCharge;
	}
	
//	public LazyEvaluateLinkedSpectrum(Peptide p1, Peptide p2, int linkedCharge){
//		this.p1 = p1;
//		this.p2 = p2;
//		this.linkedCharge = linkedCharge;
//		this.parentMass = ((p1.getParentmass()-p1.getCharge()*Mass.PROTON_MASS)+linkedCharge*Mass.PROTON_MASS)/linkedCharge;
//		this.charge = linkedCharge;
//		this.peptide = p1 + " & " + p2;
//	}
	
	
	public LazyEvaluateLinkedSpectrum(Peptide p1, Peptide p2, int linkedCharge, double linkerMass){
		this.p1 = p1;
		this.p2 = p2;
		this.linkedCharge = linkedCharge;
		this.linkerMass = linkerMass;
		this.parentMass = ((p1.getParentmass()-p1.getCharge()*Mass.PROTON_MASS)+linkedCharge*Mass.PROTON_MASS)/linkedCharge;
		this.charge = linkedCharge;
		this.peptide = p1 + " & " + p2;
	}
	
	@Override
	public List<Peak> getPeak(){
		if(p1 == null || p2 == null){
			TheoreticalSpectrum t = new TheoreticalSpectrum(this.p, linkedCharge);
//			System.out.println("generated peptide: " + t.p);
//			for(int i = 0; i < t.getPeak().size(); i++){
//				System.out.println(t.getPeak().get(i));
//			}
			return t.getPeak();
		}else{
			//System.out.println("generating peptide: " + this.peptide);
			//System.out.println(p1 + "\t" + p2);
			TheoreticalSpectrum t = new TheoreticalSpectrum(this.p1, this.p2, (short)linkedCharge, true, linkerMass);
//			System.out.println("generated peptide: " + t.p);
//			for(int i = 0; i < t.getPeak().size(); i++){
//				System.out.println(t.getPeak().get(i));
//			}
			return t.getPeak();
		}
	}
	
	public List<Peak> getPeaks(){
		return getPeak();
	}
	
	public void createSpectrum(){
		TheoreticalSpectrum t;
		if(p1 == null || p2 == null){
			t = new TheoreticalSpectrum(this.p, linkedCharge);
		}else{
			t = new TheoreticalSpectrum(this.p1, this.p2, (short)linkedCharge, true, 
					this.linkerMass);
		}
		this.p = t.p;
		this.parentMass = t.parentMass;
		this.charge = t.charge;
		this.peptide = t.peptide;
		this.setPeaks(t.getPeak());
	}
	
	public Spectrum createArraySpectrum(){
		Spectrum t = null;
		if(p1 == null || p2 == null){
			t = LinkedTheoSpectrumFactory.getLinkedTheoSpectrum(p.getPeptide(), charge, 
					p.getPos(), p.getPtmmasses(), p.getLinkedPos(), 
					LinkedTheoSpectrumFactory.linkedTypeMap, LinkedTheoSpectrumFactory.linkedIonMap);
			//System.out.println("peptide is: " + p);
			//t.peptide = p.getPeptide();
			t.peptide = p.toString();
		}else{
			int charge = this.linkedCharge;
			t = LinkedTheoSpectrumFactory.getLinkedTheoSpectrum(p1.getPeptide(), p2.getPeptide(), charge, charge, 
				 p1.getPos(), p2.getPos(), p1.getPtmmasses(), p2.getPtmmasses(), p1.getLinkedPos(), p2.getLinkedPos(), 
				 LinkedTheoSpectrumFactory.linkedTypeMap, LinkedTheoSpectrumFactory.linkedIonMap);
			t.peptide = p1.getPeptide() + " & " + p2.getPeptide();
		}
		//this.p = t.p;
		t.parentMass = p1.getParentmass();
		t.charge = p1.getCharge();
		return t;
	}

}
