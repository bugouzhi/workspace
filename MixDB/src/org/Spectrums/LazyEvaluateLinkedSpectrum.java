package org.Spectrums;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

//lazy-evaluated version of theoretical spectrum
//only generate the spectrum upon requested
//this is so to handle large spectrum library
public class LazyEvaluateLinkedSpectrum extends TheoreticalSpectrum{
	private int linkedCharge;
	Peptide p1;
	Peptide p2;
	public int getLinkedCharge() {
		return linkedCharge;
	}
	
	public void setLinkedCharge(int linkedCharge) {
		this.linkedCharge = linkedCharge;
	}
	
	
	//think about what is proper way to do this? is it to initialize everything but not created the peasks?
	public LazyEvaluateLinkedSpectrum(Peptide p, int linkedCharge){
		this.linkedCharge = linkedCharge;
		this.peptide = p.getPeptide()+"."+p.getCharge();
		this.p = p;
		this.charge = linkedCharge;
		this.parentMass = ((p.getParentmass()-p.getCharge()*Mass.PROTON_MASS)+linkedCharge*Mass.PROTON_MASS)/linkedCharge;
	}
	
	public LazyEvaluateLinkedSpectrum(Peptide p1, Peptide p2, int linkedCharge){
		this.p1 = p1;
		this.p2 = p2;
		this.linkedCharge = linkedCharge;
		this.parentMass = ((p1.getParentmass()-p1.getCharge()*Mass.PROTON_MASS)+linkedCharge*Mass.PROTON_MASS)/linkedCharge;
		this.charge = linkedCharge;
		this.peptide = p1 + " & " + p2;
	}
	
	@Override
	public List<Peak> getPeak(){
		if(p1 == null || p2 == null){
			TheoreticalSpectrum t = new TheoreticalSpectrum(this.p, linkedCharge);
			return t.getPeak();
		}else{
			TheoreticalSpectrum t = new TheoreticalSpectrum(this.p1, this.p2, (short)linkedCharge, true);
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
			t = new TheoreticalSpectrum(this.p1, this.p2, (short)linkedCharge, true);
		}
		this.p = t.p;
		this.parentMass = t.parentMass;
		this.charge = t.charge;
		this.peptide = t.peptide;
		this.setPeaks(t.getPeak());
	}

}
