package org.Spectrums;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import Utils.FileIOUtils;

/**
 * This class provide general methods for analyze peptides
 * @author jian wang
 *
 */

public class PeptideMassAnalysis {
	//private static MassCalc mcalc = initialize();
	/**
	 * Compute molecular weights for a peptide
	 * @param peptide
	 * @return mass of the peptide
	 */
	
	public static double computeMolecularMass(String pep){
		double mass = 0.0;
		//mass = mcalc.getMass(ProteinTools.createProtein(pep));			
		//mass = mcalc.getMolecularWeight(ProteinTools.createProtein(pep));
		String seq = pep;
		for(int i = 0; i < seq.length(); i++){
			mass += Mass.getAAMass(seq.charAt(i));
		}
		mass+= Mass.WATER;
		return mass;
	}
	
	public static double computeMbyZ(String pep, int charge){
		String[] tokens = pep.split("\\.");
		//System.out.println("peptide original is: " + pep);
		String fixed = tokens[0].replaceAll("\\d|\\+", "");
		//System.out.println("peptide fixed is: " + fixed);
		double mass = computeMolecularMass(fixed);
		return (mass + charge*1.007282675)/charge;
	}
	
	public static double computeMbyZ(Peptide pep, int charge){
		double mass = 0;
		if(pep.getPeptide().equals("Z")){
			mass = Mass.WATER;
		}else{
			String seq = pep.getPeptide();
			for(int i = 0; i < seq.length(); i++){
				mass += Mass.getAAMass(seq.charAt(i));
			}
			mass+= Mass.WATER;
		}
		if(pep.hasPTMs()){
			for(int i = 0; i < pep.getPtmmasses().length; i++){
				mass+= pep.getPtmmasses()[i];
			}
		}
		return (mass + charge*1.007282675)/charge;
	}
	
	//given a mass spectrum library and its annotation
	//compute the error between theoretical values and
	//the actual observed values
	public static void parentMassError(String spectrumFile, String annotationFile){
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
		System.out.println("finish loading");
		lib1.annoateSpectrumFromInspectFile(annotationFile);
		System.out.println("finish annotating");
		Iterator<Spectrum> it = lib1.iterator();
		Spectrum s;
		Peptide p;
		while(it.hasNext()){
			s = it.next();
			if(Peptide.isValidPeptide(s.peptide) && s.peptide.contains("+")){
				p = new Peptide(s.peptide);
				System.out.println(s.spectrumName + " with mass error " +
						(s.parentMass-computeMbyZ(p, p.getCharge())));
			}
		}
		
	}
	
	public static void testMassComputation(){
		String pep = "YIASLKQR";
		int c = 2;
		System.out.println("peptide is: " + pep + " weights: " + computeMolecularMass(pep) + " Dalton");
		System.out.println("peptide is: " + pep + " with charge:  " + c + " weights " + computeMbyZ(pep, c) + " Dalton");

	}
	
	public static void testparentMassError(){
		String spectrumFile = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
		String annotation =".\\mixture_linked\\trps\\result.txt";
		parentMassError(spectrumFile, annotation);
	}
	
	public static void testLinkedMass(){
		String pep1 = "TPAKEIDR";
		String pep2 = "NPYKFISR	";
		LinkedPeptide p = new LinkedPeptide(pep1 + "--" + pep2, 3);
		System.out.println("precurosr mass is: " + p.getParentmass());
	}
	
	public static void main(String[] args){
			testMassComputation();
			//testLinkedMass();
			//testparentMassError();
	}
	
	

}
