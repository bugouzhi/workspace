package org.Spectrums;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * This class help create modified version of a spectrum by
 * allowing for a mass shift in the peptide
 * @author Jian
 *
 */
public class ModifiedSpectrum {
	Spectrum s;
	public ModifiedSpectrum(Spectrum s){
		this.s = s;
	}
	
	public List<Spectrum> shiftSpectrum(char c, double mass){
		List<Spectrum> modified = new ArrayList<Spectrum>();
		String pep = s.peptide;
		for(int i = 0; i < pep.length(); i++){
			if(pep.charAt(i) == c){
				Spectrum mod = shiftSpectrum(i+1, mass);
				if(mod != null) modified.add(mod);
			}
		}
		return modified;
	}
	
	public List<Spectrum> shiftSpectrum(List<Integer> pos, double mass){
		List<Spectrum> modified = new ArrayList<Spectrum>();
		String pep = s.peptide;
		for(int i = 0; i < pos.size(); i++){
			Spectrum mod = shiftSpectrum(pos.get(i)+1, mass);
			//System.out.println("mod spec: " + mod.peptide);
			if(mod != null) modified.add(mod);
		}
		return modified;
	}
	
	public Spectrum shiftSpectrum(int pos, double mass){
		TheoreticalSpectrum t = new TheoreticalSpectrum(s.peptide+"."+s.charge);
		Spectrum modCopy = new Spectrum(s);
		SimpleMatchingGraph g = t.getMatchGraph(modCopy, 0.05);
		Collection<Peak> observed = g.getVerticeWithEdges(SimpleMatchingGraph.Observed, 1);
		int shiftCount = 0;
		//System.out.println("pos is: " + pos);
		for(Iterator<Peak> it = observed.iterator(); it.hasNext();){
			Peak p = it.next();
			LabelledPeak annot = (LabelledPeak)g.getNeighbors(p).iterator().next();
			//System.out.println("peak:\t" + p + "\t" + annot);
			if(annot.isPrefixPeak() && annot.getPos() >= pos){
				//System.out.println("Shfit: " + p + "\t" + annot);
				p.setMoz(p.getMass()+mass/annot.getCharge());
				shiftCount++;
			}
			if(annot.isSuffixPeak() && annot.getPos() > s.peptide.length() - pos){
				//System.out.println("Shfit: " + p + "\t" + annot);
				p.setMoz(p.getMass()+mass/annot.getCharge());
				shiftCount++;
			}
			
		}
	//	System.out.println("Total number of shift peaks: " + shiftCount + " total: " + this.s.getPeak().size());
		Peptide p = new Peptide(s.peptide, s.charge);
		p.insertPTM(pos, mass);
		//System.out.println("mod pep " + p);
		modCopy.peptide = p.toString();
		modCopy.parentMass = p.getParentmass();
		if(checkModSpec(this.s, modCopy)){
			return modCopy;
		}else{
			return null;
		}
	}
	
	public static void addModToLibrary(List<Spectrum> specList){
		int size = specList.size();
		for(int i = 0; i < size; i++){
			Spectrum s = specList.get(i);
			if(s.peptide.contains("+")){
				continue;
			}
			ModifiedSpectrum mod = new ModifiedSpectrum(s);
			//Spectrum s2 = mod.shiftSpectrum(0, 42.011);
			//if(s2!=null){specList.add(s2);};
			if(s.peptide.contains("M")){
				specList.addAll(mod.shiftSpectrum('M', 15.9999));				
			}
		}
		System.out.println("Library size: " + size + "\tmoded\t" + specList.size());
		for(int i = 0; i < specList.size(); i++){
			
		}
	}
	
	/**
	 * currently support single PTM only (possibly extend to multi-ptms)
	 * @param specList
	 * @param ptms
	 */
	public static void addModToLibrary(List<Spectrum> specList, List<PTM[]> ptms){
		int size = specList.size();
		for(int i = 0; i < ptms.size(); i++){
			PTM[] ptmCom = ptms.get(i);
			PTM ptm = ptmCom[0];
			for(int j = 0; j < size; j++){
				Spectrum s = specList.get(j);
				if(s.peptide.contains("+")){
					continue;
				}
				ModifiedSpectrum mod = new ModifiedSpectrum(s);
				List<Integer> posList = ptm.getPTMPositions(s.peptide);
				specList.addAll(mod.shiftSpectrum(posList, ptm.getPtmMass()));
			}
		}
		System.out.println("Library size: " + size + "\tmoded\t" + specList.size());
	}
	
	private boolean checkModSpec(Spectrum s, Spectrum mod){
		//System.out.println(s.peptide + "\t" + mod.peptide + "\t" + s.charge +  " Sim: " + s.cosine(mod, 0.05));
		if(s.cosine(mod, 0.05) > 0.85){
			return false;
		}else{
			return true;
		}
	}
	
	public static void testModifiedSpectrum(){
		String libFile = "..//mixture_linked//ACG_swathdevelopment_UPSEcoli_REP234_IDA_plusDecoy2.mgf";
		SpectrumLib lib = new SpectrumLib(libFile, "MGF");
		addModToLibrary(lib.getAllSpectrums());
	}
	
	public static void main(String[] args){
		testModifiedSpectrum();
	}
	
}
