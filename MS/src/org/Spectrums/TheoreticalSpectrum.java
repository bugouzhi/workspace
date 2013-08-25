package org.Spectrums;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Arrays;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.util.TreeSet;
import java.util.HashMap;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import Utils.FileIOUtils;



/**
 * Theoretical Spectrum of a peptide
 * @author jian wang
 *
 */
public class TheoreticalSpectrum extends Spectrum{
	public static int MAX_PROFILE_SPAN = 0;
	protected static String[] prefixIons = {"b", "b-H20", "b-NH3", "b(iso)"};//, "a", "a-H20", "a-NH3"};
	protected static  String[] suffixIons = {"y", "y-H20", "y-NH3", "y(iso)"};
	public static boolean deconvolutedMode=false;
	//protected static String[] prefixIons = Mass.standardPrefixes;
	//protected static  String[] suffixIons = Mass.standardSuffixes;
	protected static double MS2Tolerance = 0.5;
	protected static int minMatchingPeaks = 10;
	
	
	public  Peptide p;
	
	
	protected TheoreticalSpectrum(){
		
	}
	
	public TheoreticalSpectrum(String peptide){
		this(new Peptide(peptide));
	}
	
	public TheoreticalSpectrum(String peptide, String[] prefixIons, String[] suffixIons){
		this(new Peptide(peptide), prefixIons, suffixIons);
	}
	
	public TheoreticalSpectrum(Peptide p){
		super();
		Vector<Peak> theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
				p.getPos(), p.getPtmmasses(), p.getCharge());
//		Vector theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
//				null, null, p.getCharge());
//		this.shiftPeaks(theoPeaks, p.getPos(), p.getPtmmasses());
		this.setPeaks(theoPeaks);
		this.peptide = p.getPeptide()+"."+p.getCharge();
		this.parentMass= PeptideMassAnalysis.computeMbyZ(p, p.getCharge());
		this.charge = p.getCharge();
		this.setPeptide(theoPeaks, p);
		this.p = p;
	
	}
	
	public static void addIsotopicPeaks(TheoreticalSpectrum th, int numC13){
		List<Peak> pList = th.getPeak();
		List<Peak> isoList = new ArrayList<Peak>(pList.size()*numC13);
		double massOffset = Mass.C13 - Mass.C12;
		for(int i = 0; i < pList.size(); i++){
			LabelledPeak curr = (LabelledPeak)pList.get(i);
			for(int c13 = 1; c13 <= numC13; c13++){
				Peak isoPeak = new LabelledPeak(curr.getMass() + (massOffset*c13)/(double)curr.getCharge(), 
						curr.getIntensity(), curr.getType()+"(iso)"+c13, curr.getPos(), curr.getCharge());
				isoList.add(isoPeak);
			}
		}
		System.out.println("we start with peaks " + pList.size());
		pList.addAll(isoList);
		Collections.sort(pList, PeakMassComparator.comparator);
		System.out.println("we end with peaks: " + th.getPeak().size());
		th.setPeaks(pList);
	}
	
	public TheoreticalSpectrum(Peptide p, String[] prefixIons, String[] suffixIons){
		super();
		p = new Peptide(p);
		this.peptide = p.getPeptide()+"."+p.getCharge();
		if(p.getCharge() == 4){
			//p.setCharge((short)3);
		}
		Vector<Peak> theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
				p.getPos(), p.getPtmmasses(), p.getCharge());
//		Vector theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
//				null, null, p.getCharge());
//		this.shiftPeaks(theoPeaks, p.getPos(), p.getPtmmasses());
		this.setPeaks(theoPeaks);
		//this.peptide = p.getPeptide()+"."+p.getCharge();
		this.parentMass= PeptideMassAnalysis.computeMbyZ(p, p.getCharge());
		this.charge = p.getCharge();
		this.setPeptide(theoPeaks, p);
		this.p = p;
	}
	
//	public TheoreticalSpectrum(String peptide1, String peptide2){
//		super();
//		Peptide p1 = new Peptide(peptide1);
//		Peptide p2 = new Peptide(peptide2);
//		this.peptide = p1.getPeptide() + "&" + p2.getPeptide();
//		short charge = (short)(p1.getCharge() + p2.getCharge());
//		double[] ptmmass1 = new double[] {PeptideMassAnalysis.computeMbyZ(p2.getPeptide(), 1)
//				+ Mass.DSPLINKER_MASS};
//		double[] ptmmass2 = new double[] {PeptideMassAnalysis.computeMbyZ(p1.getPeptide(), 1)
//				+ Mass.DSPLINKER_MASS};
//		//ptmmass1[0] = 0;
//		//ptmmass2[0] = 0;
//
//		Vector<Peak> theoPeaks = this.generatePeaks(p1.getPeptide(), prefixIons, suffixIons,
//				p1.getPos(), ptmmass1, charge);
//		Vector<Peak> theoPeaks2 = this.generatePeaks(p2.getPeptide(), prefixIons, suffixIons,
//				p2.getPos(), ptmmass2, charge);
//		setPeptide(theoPeaks, p1);
//		setPeptide(theoPeaks2, p2);
//		theoPeaks.addAll(theoPeaks2);
//		Collections.sort(theoPeaks, new PeakMassComparator());
//		System.out.println("total peaks created: " + theoPeaks.size());		
//		this.setPeaks(theoPeaks);
//		this.parentMass= (ptmmass1[0] + ptmmass2[0]+Mass.DSPLINKER_MASS + charge*Mass.PROTON_MASS)/charge;
//		this.charge = charge;
//	}
	
	//for mixture spectrum
	
	
	public TheoreticalSpectrum(String peptide1, String peptide2){
		super();
		Peptide p1 = new Peptide(peptide1);
		Peptide p2 = new Peptide(peptide2);
		//this.peptide = p1.getPeptide() + " & " + p2.getPeptide();
		this.peptide = p1.toString() + " & " + p2.toString();
		if(p1.getCharge()== 4){
			//p1.setCharge((short)3);
		}
		if(p2.getCharge() == 4){
			//p2.setCharge((short)3);
		}
		short charge = (short)(p1.getCharge() + p2.getCharge());
		double[] ptmmass1 = new double[] {0};
		double[] ptmmass2 = new double[] {0};
		//ptmmass1[0] = 0;
		//ptmmass2[0] = 0;

		Vector<Peak> theoPeaks = this.generatePeaks(p1.getPeptide(), prefixIons, suffixIons,
				p1.getPos(), ptmmass1, p1.getCharge());
		Vector<Peak> theoPeaks2 = this.generatePeaks(p2.getPeptide(), prefixIons, suffixIons,
				p2.getPos(), ptmmass2, p2.getCharge());
		Vector<Peak> mixturePeaks = new Vector<Peak>();
		setPeptide(theoPeaks, p1);
		setPeptide(theoPeaks2, p2);
		
		for(int i = 0; i < theoPeaks.size(); i++){
			MixturePeak m = new MixturePeak((LabelledPeak)theoPeaks.get(i), 0);
			m.setParent(this);
			mixturePeaks.add(m);
		}
		
		for(int i = 0; i < theoPeaks2.size(); i++){
			MixturePeak m = new MixturePeak((LabelledPeak)theoPeaks2.get(i), 1);
			m.setParent(this);
			mixturePeaks.add(m);			
		}
		Collections.sort(mixturePeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());		
		this.setPeaks(mixturePeaks);
		this.peptide = peptide1 + " & " + peptide2;
		this.parentMass= (ptmmass1[0] + ptmmass2[0])/2;
		this.charge = charge;
	}
	
	//for linked
	public TheoreticalSpectrum(Spectrum s1, Spectrum s2){
		super();
		String peptide1 = s1.peptide, peptide2 = s2.peptide;
		Peptide p1 = new Peptide(peptide1);
		Peptide p2 = new Peptide(peptide2);

		this.peptide = p1.getPeptide() + "&" + p2.getPeptide();
		short charge = (short)(p1.getCharge() + p2.getCharge());
		double[] ptmmass1 = new double[] {PeptideMassAnalysis.computeMbyZ(p2.getPeptide(), 1)
				+ Mass.DSPLINKER_MASS - 146.0}; //we discount an extra hydrogen here since it is include in weight of 2nd peptide
		double[] ptmmass2 = new double[] {PeptideMassAnalysis.computeMbyZ(p1.getPeptide(), 1)
				+ Mass.DSPLINKER_MASS - 146.0};
		Vector<Peak> theoPeaks1 = this.generatePeaks(p1.getPeptide(), prefixIons, suffixIons,
				p1.getPos(), new double[]{145}, p1.getCharge());
		Vector<Peak> theoPeaks2 = this.generatePeaks(p2.getPeptide(), prefixIons, suffixIons,
				p2.getPos(), new double[]{145}, p2.getCharge());
		
		TheoreticalSpectrum th1 = new TheoreticalSpectrum(s1.peptide);
		TheoreticalSpectrum th2 = new TheoreticalSpectrum(s2.peptide);

		List matchedPeaks1 = th1.getPeak();//annotatePeaks(th1, s1);
		List matchedPeaks2 = th2.getPeak();//annotatePeaks(th2, s2);
		Vector<Peak> theoPeaks = new Vector();

		
//		filteredList1[0] = filteredList1[1];//theoPeaks1; //
//		filteredList2[0] = filteredList2[1];// theoPeaks2; //
		this.shiftPeaks(matchedPeaks1, p1.getPos(), ptmmass1, p1.getPeptide().length());
		this.shiftPeaks(matchedPeaks2, p2.getPos(), ptmmass2, p2.getPeptide().length());

		this.addChargedPeaks(matchedPeaks1, p1.getPos(), p1.getCharge()+p2.getCharge());
		this.addChargedPeaks(matchedPeaks2, p2.getPos(), p1.getCharge()+p2.getCharge());
		
		setPeptide(matchedPeaks1, p1);
		setPeptide(matchedPeaks2, p2);
		
		theoPeaks.addAll(matchedPeaks1);
		theoPeaks.addAll(matchedPeaks2);

		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		this.parentMass= (s1.parentMass + s2.parentMass+Mass.DSPLINKER_MASS + charge*Mass.PROTON_MASS)/charge;
		this.charge = charge;
	}
	
	//for linked peptides
	public TheoreticalSpectrum(Peptide p1, int linkedCharge){
		//this(p1, linkedCharge, new String[]{"b", "b-H20", "b-NH3"}, new String[]{"y", "y-H20", "y-NH3"});
		this(p1, linkedCharge, TheoreticalSpectrum.prefixIons, TheoreticalSpectrum.suffixIons);
		// this(p1, linkedCharge, new String[]{"b"}, new String[]{"y"});
	}
	
	public TheoreticalSpectrum(Peptide p1, int linkedCharge, String[] prefix, String[] suffix){
		this.peptide = p1.getPeptide()+"." + linkedCharge;
		p1.setCharge((short)linkedCharge);
		List theoPeaks = this.generatePeaks(p1.getPeptide(), prefix, suffix, p1.getPos(), p1.getPtmmasses(), p1.getCharge());
		//p1.setCharge((short)1);
		//this.shiftPeaks(theoPeaks, p1.getPos(), p1.getPtmmasses(), p1.getPeptide().length());
		//this.addChargedPeaks(theoPeaks, p1.getPos(), linkedCharge);
		this.setLinkedPeptide(theoPeaks, p1, linkedCharge);
		p1.setCharge((short)linkedCharge);
		this.setPeptide(theoPeaks, p1);	
		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println(this.peptide + " total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		this.parentMass = p1.getParentmass(); 
		//this.charge = p1.getCharge();
		this.charge = linkedCharge;
		this.p = p1;
		this.peptide = p1.getPeptide() + "." + p1.getCharge();
	}
	
	//for linked peptide
	public TheoreticalSpectrum(Peptide p1, Peptide p2, short charge){
		super();
		this.peptide = p1.getPeptide() + "&" + p2.getPeptide();
		double[] ptmmass1 = new double[] {PeptideMassAnalysis.computeMbyZ(p2.getPeptide(), 1)
				+ Mass.DSPLINKER_MASS -146 }; //we discount an extra hydrogen here since it is include in weight of 2nd peptide
		double[] ptmmass2 = new double[] {PeptideMassAnalysis.computeMbyZ(p1.getPeptide(), 1)
				+ Mass.DSPLINKER_MASS -146 };
		
		TheoreticalSpectrum th1 = new TheoreticalSpectrum(p1);
		TheoreticalSpectrum th2 = new TheoreticalSpectrum(p2);

		List matchedPeaks1 = th1.getPeak();//annotatePeaks(th1, s1);
		List matchedPeaks2 = th2.getPeak();//annotatePeaks(th2, s2);
		Vector<Peak> theoPeaks = new Vector();

		this.shiftPeaks(matchedPeaks1, p1.getPos(), ptmmass1, p1.getPeptide().length());
		this.shiftPeaks(matchedPeaks2, p2.getPos(), ptmmass2, p2.getPeptide().length());

		this.addChargedPeaks(matchedPeaks1, p1.getPos(), charge);
		this.addChargedPeaks(matchedPeaks2, p2.getPos(), charge);
		
		p1.setCharge(charge);
		p2.setCharge(charge);
		setPeptide(matchedPeaks1, p1);
		setPeptide(matchedPeaks2, p2);
		
		theoPeaks.addAll(matchedPeaks1);
		theoPeaks.addAll(matchedPeaks2);

		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		this.parentMass= (p1.getParentmass() + p2.getParentmass() + Mass.DSPLINKER_MASS + charge*Mass.PROTON_MASS)/charge;
		this.charge = charge;
	}
	
	public TheoreticalSpectrum(Peptide p1, Peptide p2, short charge, boolean dummy){
		super();
		this.peptide = p1 + " & " + p2;
		//System.out.println("p1: " + p1 + "\tp2: " + p2);
		this.p = new LinkedPeptide(p1, p2, charge);
//		System.out.println("Peptide pairs are: " + p1 + "\t" + p1.getCharge() + "\t" + p1.getParentmass() + "\t" + p1.getLinkedPos() + "\t"
//				+ p2 + "\t" + p2.getCharge() + "\t" + p2.getParentmass() + "\t" + p2.getLinkedPos());
		TheoreticalSpectrum th1 = new TheoreticalSpectrum(p1, charge, this.prefixIons, this.suffixIons);
		TheoreticalSpectrum th2 = new TheoreticalSpectrum(p2, charge, this.prefixIons, this.suffixIons);
//		TheoreticalSpectrum th1 = new TheoreticalSpectrum(p1, charge, this.prefixIons, this.suffixIons);
//		TheoreticalSpectrum th2 = new TheoreticalSpectrum(p2, charge, this.prefixIons, this.suffixIons);


		List matchedPeaks1 = th1.getPeak();//annotatePeaks(th1, s1);
		List matchedPeaks2 = th2.getPeak();//annotatePeaks(th2, s2);
		List mixturePeaks1 = new ArrayList();//annotatePeaks(th1, s1);
		List mixturePeaks2 = new ArrayList();//annotatePeaks(th2, s2);

		Vector<Peak> theoPeaks = new Vector();
		
		for(int i = 0; i < matchedPeaks1.size(); i++){
			MixturePeak p = new MixturePeak((LabelledPeak)matchedPeaks1.get(i), 0);
			p.setParent(this);
			mixturePeaks1.add(p);
		}
		this.setPeptide(mixturePeaks1, p1);
		
		for(int i = 0; i < matchedPeaks2.size(); i++){
			MixturePeak p = new MixturePeak((LabelledPeak)matchedPeaks2.get(i), 1);
			p.setParent(this);
			mixturePeaks2.add(p);
		}
		this.setPeptide(mixturePeaks2, p2);
		
//		if(charge <= 4){
//			p1.setCharge((short)2);
//			p2.setCharge((short)2);
//		}else{
//			p1.setCharge((short)3);
//			p2.setCharge((short)3);
//		}
		//setLinkedPeptide(matchedPeaks1, p1, charge);
		//setLinkedPeptide(matchedPeaks2, p2, charge);
		
		theoPeaks.addAll(mixturePeaks1);
		theoPeaks.addAll(mixturePeaks2);

		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		//LinkedPeptide lp = new LinkedPeptide(p1.getPeptide() + "--" + p2.getPeptide(), charge);
		this.parentMass= this.p.getParentmass();
		this.charge = charge;
//		for(int i = 0; i < theoPeaks.size(); i ++){
//			System.out.println(theoPeaks.get(i));
//		}
	}
	
	public TheoreticalSpectrum(LinkedPeptide lp, short charge, boolean dummy){
		super();
		this.peptide = lp.peptides[0] + " & " + lp.peptides[1];
		
		TheoreticalSpectrum th1 = new TheoreticalSpectrum(lp.peptides[0], charge, this.prefixIons, this.suffixIons);
		TheoreticalSpectrum th2 = new TheoreticalSpectrum(lp.peptides[1], charge, this.prefixIons, this.suffixIons);
//		TheoreticalSpectrum th1 = new TheoreticalSpectrum(p1, charge, this.prefixIons, this.suffixIons);
//		TheoreticalSpectrum th2 = new TheoreticalSpectrum(p2, charge, this.prefixIons, this.suffixIons);

		List matchedPeaks1 = th1.getPeak();//annotatePeaks(th1, s1);
		List matchedPeaks2 = th2.getPeak();//annotatePeaks(th2, s2);
		List mixturePeaks1 = new ArrayList();//annotatePeaks(th1, s1);
		List mixturePeaks2 = new ArrayList();//annotatePeaks(th2, s2);

		Vector<Peak> theoPeaks = new Vector();
		
		for(int i = 0; i < matchedPeaks1.size(); i++){
			MixturePeak p = new MixturePeak((LabelledPeak)matchedPeaks1.get(i), 0);
			p.setParent(this);
			mixturePeaks1.add(p);
		}
		this.setPeptide(mixturePeaks1, lp.peptides[0]);
		
		for(int i = 0; i < matchedPeaks2.size(); i++){
			MixturePeak p = new MixturePeak((LabelledPeak)matchedPeaks2.get(i), 1);
			p.setParent(this);
			mixturePeaks2.add(p);
		}
		this.setPeptide(mixturePeaks2, lp.peptides[1]);
		
//		if(charge <= 4){
//			p1.setCharge((short)2);
//			p2.setCharge((short)2);
//		}else{
//			p1.setCharge((short)3);
//			p2.setCharge((short)3);
//		}
		//setLinkedPeptide(matchedPeaks1, p1, charge);
		//setLinkedPeptide(matchedPeaks2, p2, charge);
		
		theoPeaks.addAll(mixturePeaks1);
		theoPeaks.addAll(mixturePeaks2);

		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		this.p = lp;
		this.parentMass= lp.getParentmass();
		this.charge = charge;
//		for(int i = 0; i < theoPeaks.size(); i ++){
//			System.out.println(theoPeaks.get(i));
//		}
	}
	
	
	public TheoreticalSpectrum(TheoreticalSpectrum t1, TheoreticalSpectrum t2, short charge, boolean dummy){
		super();
		this.peptide = t1.getPeptide() + " & " + t2.getPeptide();
		
		TheoreticalSpectrum th1 = t1;
		TheoreticalSpectrum th2 = t2;
//		TheoreticalSpectrum th1 = new TheoreticalSpectrum(p1, charge, this.prefixIons, this.suffixIons);
//		TheoreticalSpectrum th2 = new TheoreticalSpectrum(p2, charge, this.prefixIons, this.suffixIons);


		List matchedPeaks1 = th1.getPeak();//annotatePeaks(th1, s1);
		List matchedPeaks2 = th2.getPeak();//annotatePeaks(th2, s2);
		Vector<Peak> theoPeaks = new Vector();

		//setLinkedPeptide(matchedPeaks1, p1, charge);
		//setLinkedPeptide(matchedPeaks2, p2, charge);
		
		theoPeaks.addAll(matchedPeaks1);
		theoPeaks.addAll(matchedPeaks2);

		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		this.parentMass= t1.parentMass;
		this.charge = charge;
	}
	public static TheoreticalSpectrum getLinkedTheoreticalSpectrum(String pep1, String pep2, short charge, int linkedPosition1, int linkedPosition2){
		Peptide p1 = new Peptide(pep1, 1);
		Peptide p2 = new Peptide(pep2, 1);
		double mass = (p1.getParentmass() + p2.getParentmass()
				+ Mass.DSSLINKER_MASS 
				+ Mass.PROTON_MASS*(charge-2))/charge;	
		double massShift1 = p1.getParentmass();
		double massShift2 = p2.getParentmass();
 		p1.insertPTM(p1.getPeptide().indexOf('K')+1, massShift2+Mass.DSSLINKER_MASS-Mass.PROTON_MASS);
		p2.insertPTM(p2.getPeptide().indexOf('K')+1, massShift1+Mass.DSSLINKER_MASS-Mass.PROTON_MASS);
		
		return new TheoreticalSpectrum(p1, p2, charge, true);
	}
	
	protected void setPeptide(List peaks, Peptide p){
		for(int i = 0; i < peaks.size(); i++){
			((LabelledPeak)peaks.get(i)).setPep(p);
		}
	}
	
	protected void setLinkedPeptide(List peaks, Peptide pep, int linkedCharge){
		Peptide linked = new Peptide(pep);
		linked.setCharge((short)linkedCharge);
		for(int i = 0; i < peaks.size(); i++){
			LabelledPeak lp = ((LabelledPeak)peaks.get(i));
			if(this.isLinkedPeak(pep, lp)){
				lp.setPep(linked);
			}else{
				lp.setPep(pep);
			}
		}
	}
	
	public String getPeptide(){
		return this.peptide;
		//return this.p.getPeptide() + "." + this.p.getCharge();
	}
	
	public List<Peak> annotatePeaks(TheoreticalSpectrum th, Spectrum s){
		List[] matchedPeaks = th.matchSpectrum(s, TheoreticalSpectrum.MS2Tolerance);
		List<Peak> annotatedPeaks = new ArrayList();
		LabelledPeak p1, p2;
		for(int i = 1; i < matchedPeaks[0].size(); i++){
			p1 = ((LabelledPeak)matchedPeaks[0].get(i-1));
			p2 = ((LabelledPeak)matchedPeaks[0].get(i));
			//if(!p1.equals(p2)){
				annotatedPeaks.add(p2.transferLabel((Peak)(matchedPeaks[1]).get(i)));
			//}
		}
		return annotatedPeaks;
	}
	
	public Vector<Peak> generatePeaks(String peptide, String[] prefixIons, String[] suffixIons){
		return generatePeaks(peptide, prefixIons, suffixIons, null, null, (short)1);
	}
	
	public Vector<Peak> generatePeaks(String peptide, String[] prefixIons, String[] suffixIons, 
			int[] pos, double[] ptmMass, short charge){
		if(charge==4){
			//charge=3;
		}
		return generatePeaks(peptide, prefixIons, suffixIons, pos, ptmMass, 1, charge);
		
	}
	
	public Vector<Peak> generatePeaks(String peptide, String[] prefixIons, String[] suffixIons,
			int[] pos, double[] ptmMass, int minCharge, int maxCharge){
		double[][] base = computeBaseMass(peptide, pos, ptmMass);
		Vector<Peak> theoPeaks = new Vector();
		//generating prefix ionsg
		double[] ions = new double[base[0].length];
		double[] chargedIons = new double[base[0].length];
		for(int i = 0; i < prefixIons.length; i++){
			this.addIonMod(base[0], ions, prefixIons[i], 0);
			for(short c = (short)minCharge; c <= maxCharge; c++){
				this.addCharge(ions, chargedIons, c, 0);
				createPeaks(theoPeaks, chargedIons, prefixIons[i], c);
			}
		}
		
//		//generating suffix ions
		for(int i = 0; i < suffixIons.length; i++){
			this.addIonMod(base[1], ions, suffixIons[i], 0);
			for(short c = (short)minCharge; c <= maxCharge; c++){
				this.addCharge(ions, chargedIons, c, 0);
				createPeaks(theoPeaks, chargedIons, suffixIons[i], c);
			}
		}
		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total ions created is: " + theoPeaks.size());
		return theoPeaks;
	}
	
	public Vector<Peak> generatePeaks(double[][] base, String[] prefixIons, String[] suffixIons,
			int[] pos, double[] ptmMass, int minCharge, int maxCharge){
		
		Vector<Peak> theoPeaks = new Vector();
		//generating prefix ionsg
		double[] ions = new double[base[0].length];
		double[] chargedIons = new double[base[0].length];
		for(int i = 0; i < prefixIons.length; i++){
			this.addIonMod(base[0], ions, prefixIons[i], 0);
			for(short c = (short)minCharge; c <= maxCharge; c++){
				this.addCharge(ions, chargedIons, c, 0);
				createPeaks(theoPeaks, chargedIons, prefixIons[i], c);
			}
		}
		
//		//generating suffix ions
		for(int i = 0; i < suffixIons.length; i++){
			this.addIonMod(base[1], ions, suffixIons[i], 0);
			for(short c = (short)minCharge; c <= maxCharge; c++){
				this.addCharge(ions, chargedIons, c, 0);
				createPeaks(theoPeaks, chargedIons, suffixIons[i], c);
			}
		}
		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total ions created is: " + theoPeaks.size());
		return theoPeaks;
	}
	
	private void createPeaks(List<Peak> pList, double[] masses, String type, short charge){
		LabelledPeak p;
		for(int i = 0; i < masses.length; i++){
			//p = new LabelledPeak(masses[i], LabelledPeak.DEFAULT_INTENS, type, (short)(i+1), charge);
			p = LabelledPeakFactory.createLabelledPeak(masses[i], LabelledPeak.DEFAULT_INTENS, type, (short)(i+1), charge);
			pList.add(p);
		}
	}
	
	public double[] computeTheoreticalMasses(String peptide, String[] prefixIons, String[] suffixIons){
		double[][] base = computeBaseMass(peptide, null, null);
		double[] prefixMasses = computeIonsMass(base[0], prefixIons);
		double[] suffixMasses = computeIonsMass(base[1], suffixIons);
		double[] allMasses = ArrayUtils.mergeArray(prefixMasses, suffixMasses);
		Arrays.sort(allMasses);
		return allMasses;
	}
	
	public static double[][] computeBaseMass(String peptide, int[] pos, double[] ptmmass){
		double[] prefixMasses = new double[peptide.length()];
		double[] suffixMasses = new double [peptide.length()];
		double sum = 0.0; 
		//prefix sums
		int j = 0;
		for(int i = 0; i < peptide.length(); i++){
			if(pos != null && j < pos.length && pos[j]== i+1){
				sum += ptmmass[j];
				j++;
			}
			sum += Mass.getAAMass(peptide.charAt(i));
			prefixMasses[i] = sum;
			//System.out.println(sum);
		}
		
		
		j = 0;
		for(int i = 0; i < peptide.length()-1; i++){
			suffixMasses[i] = sum - prefixMasses[peptide.length()-i-2];

		}
		//suffixMasses[peptide.length()-1] = prefixMasses[peptide.length()-1];
		//suffixMasses[peptide.length()-1] = 0; //we zero out the unfragmented precursor masses cause it 
		prefixMasses[peptide.length()-1] = 0; //do not provide information for IDs
		return new double[][]{prefixMasses, suffixMasses};
		
	}
	
	public double[] computeIonsMass(double[] baseMass, String[] ions){
		double[] ionMasses = new double[baseMass.length*ions.length];
		int startInd = 0;
		for(int i = 0; i < ions.length; i++){
			addIonMod(baseMass, ionMasses, ions[i], startInd);
			startInd += baseMass.length;
		}
		return ionMasses;
	}
	
	public double[] computeIonsMass(double[] baseMass, double[] ionsOffSet){
		double[] ionMasses = new double[baseMass.length*ionsOffSet.length];
		int startInd = 0;
		for(int i = 0; i < ionsOffSet.length; i++){
			addIonMod(baseMass, ionMasses, ionsOffSet[i], startInd);
			startInd += baseMass.length;
		}
		return ionMasses;
	}
	
	private double[] addIonMod(double[] baseMass, double[] ionMasses, String ion, int startInd){
		double ionMass = Mass.getIonMod(ion);
		for(int i = 0; i < baseMass.length; i++){
			ionMasses[startInd+i] = baseMass[i] + ionMass;
		}
		return ionMasses;
	}
	
	private double[] addIonMod(double[] baseMass, double[] ionMasses, double ionOffSet, int startInd){
		for(int i = 0; i < baseMass.length; i++){
			ionMasses[startInd+i] = baseMass[i] + ionOffSet;
		}
		return ionMasses;
	}
	
	private double[] addCharge(double[] ionMasses, double[] chargedIonMass, int charge, int startInd){
		for(int i = 0; i < ionMasses.length; i++){
			if(this.deconvolutedMode){
				chargedIonMass[startInd+i] = (ionMasses[i] + Mass.PROTON_MASS*(charge-1));
			}else{
				chargedIonMass[startInd+i] = (ionMasses[i] + Mass.PROTON_MASS*(charge-1))/(double)charge;
			}
		}
		return ionMasses;
	}
	
	private double[] computeMultiChargedMass(){
		return new double[]{0.0};
	}
	
	public double explainedPeaks2(Spectrum theoretical, Spectrum actual, double tolerance){
		List<Peak>[] matchedPeaks = this.matchSpectrum(actual, tolerance);
		return explainedPeaks2(matchedPeaks, actual, tolerance);
	}
	
	public double explainedPeaks2(List<Peak>[] matchedPeak, Spectrum actual, double tolerance){
		return explainedPeaks2(matchedPeak, actual, tolerance, true);
	}
	
	public double explainedPeaks2(List<Peak>[] matchedPeak, Spectrum actual, double tolerance, boolean detail){
		double total = actual.sumOfPeaks();
		HashMap<Peak, Peak> presented = new HashMap();
		total = total*total;
		double explained = 0.0;
		double intense = 0.0;
		Peak prev = null;
		for(int i = 0; i < matchedPeak[1].size(); i++){
			if(!presented.containsKey(matchedPeak[1].get(i))){ //avoid double counting
				LabelledPeak lp = (LabelledPeak)matchedPeak[0].get(i);
				if(matchedPeak[1].get(i).getIntensity() > 0.15){
					intense += matchedPeak[1].get(i).getIntensity();
				}
//				if(lp.getPos() == lp.getPep().getPeptide().length()){ //subtract unfragmented precursors and its variants
//					total -= matchedPeak[1].get(i).getIntensity();
//				}
				explained += matchedPeak[1].get(i).getIntensity();
				presented.put(matchedPeak[1].get(i), matchedPeak[1].get(i));
			}
		}
		
		if(detail){
			//System.out.println("total : " + total + "\texpalined : " + explained);
			System.out.println(actual.spectrumName + "(" + actual.peptide +  " charge: " + actual.charge + ") has explained b/y ions fraction: " + 
					explained/total + " intense b/y: " + intense/total);
		}
		return explained/total;
	}

	public double[] IonsStat(List[] peaks, Spectrum s){
		return IonsStat(peaks, s, true);
	}
	public double[] IonsStat(List[] peaks, Spectrum s, boolean detail){
		if(peaks[0].size() == 0){
			return new double[]{0,0,0,0,0,0};
		}
		double[] b_y_presented = new double[6];
		List<TreeSet> presented = new ArrayList();
		TreeSet<Integer> b = new TreeSet<Integer>();
		TreeSet<Integer> y = new TreeSet<Integer>();

		LabelledPeak l;
		for(int i = 0; i < peaks[0].size(); i++){
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){	
				l = (LabelledPeak)peaks[0].get(i);
				//System.out.println("checking matched peak " + l);
				if(l.isPrefixPeak() && l.isPrimary()){ //&& l.getType().equals("b")){
					b.add(new Integer(l.getPos()));
				}
				if(l.isSuffixPeak() && l.isPrimary()){ //&& l.getType().equals("y")){
					y.add(new Integer(l.getPos()));
				}
			}
		}
		String pep = ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide();
		b_y_presented[0] = (double)b.size() / (double)pep.length();
		b_y_presented[1] = (double)y.size() / (double)pep.length();
		if(detail)
			System.out.print(s.spectrumName + "(" + ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide() + "): b/y series:" );
		Iterator<Integer> it = b.iterator();
		int prev = -1, current = 0;
		double bSeries = 0;
		double maxBSeries = 0;
		while(it.hasNext()){
			current = it.next().intValue();
			if(current - prev == 1){
				if(detail)
					System.out.print("--b"+current);
				bSeries++;
			}else{
				if(detail)
					System.out.print("  b"+current);
				bSeries = 0;
			}
			maxBSeries = bSeries > maxBSeries ? bSeries : maxBSeries;
			prev = current;
		}
		if(detail)
			System.out.println("\tfraction of total:\t" + b_y_presented[0]);
		it = y.iterator();
		prev = -1;
		current = 0;
		if(detail)
			System.out.print(s.spectrumName + "(" + ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide() + ") b/y series: " );
		double ySeries = 0;
		double maxYSeries = 0;
		while(it.hasNext()){
			current = it.next().intValue();
			if(current - prev == 1){
				if(detail)
					System.out.print("--y"+current);
				ySeries++;
			}else{
				if(detail)
					System.out.print("  y"+current);
				ySeries = 0;
			}
			maxYSeries = ySeries > maxYSeries ? ySeries : maxYSeries;
			prev = current;
		}
		if(detail)
			System.out.println( "\tfration of total:\t" + b_y_presented[1]);
		b_y_presented[2] = maxBSeries+1;
		b_y_presented[3] = maxYSeries+1;
		int[] maxConsecSeries = this.getMaxSeries(peaks, s, 1);
		b_y_presented[4] = maxConsecSeries[0];
		b_y_presented[5] = maxConsecSeries[1];
		return b_y_presented;
	}
	
	public TreeMap<Integer, Peak>[] getMatchedPeakRank(List[] peaks, Spectrum s){
		LabelledPeak l;
		TreeMap<Integer, Peak> bRanks = new TreeMap<Integer, Peak>();
		TreeMap<Integer, Peak> yRanks = new TreeMap<Integer, Peak>();
		for(int i = 0; i < peaks[0].size(); i++){
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){	
				l = (LabelledPeak)peaks[0].get(i);
				if(l.isPrefixPeak() && l.isPrimary()){ //&& l.getType().equals("b")){
					Peak p = (Peak)peaks[1].get(i);
					bRanks.put(p.getRank(), l);
				}
				if(l.isSuffixPeak() && l.isPrimary()){ //&& l.getType().equals("y")){
					Peak p = (Peak)peaks[1].get(i);
					yRanks.put(p.getRank(), l);
				}
			}
		}
		return new TreeMap[]{bRanks, yRanks};
	}
	
	public int[] getMaxRanks(int topN, List[] peaks, Spectrum s){
		TreeMap<Integer, Peak>[] matchedPeakRanks = getMatchedPeakRank(peaks, s);
		int count = 0;
		int maxRank1=0, maxRank2=0;
		System.out.print("ions ranks:\t");
		for(Iterator<Integer> it = matchedPeakRanks[0].keySet().iterator(); it.hasNext(); count++){
			Integer curr = it.next();
			LabelledPeak p = (LabelledPeak)matchedPeakRanks[0].get(curr);
			if(count <= topN) System.out.print(curr + "\t");//+ " (" + p +") ");
			if(count == topN){
				maxRank1 = curr.intValue();
			}
		
		}
		count = 0;
		for(Iterator<Integer> it = matchedPeakRanks[1].keySet().iterator(); it.hasNext(); count++){
			Integer curr = it.next();
			LabelledPeak p = (LabelledPeak)matchedPeakRanks[1].get(curr);
			if(count <= topN) System.out.print(curr + "\t");//  + " (" + p + ") ");
			if(count == topN){
				maxRank2 = curr.intValue();
			}
	
		}
		System.out.println();
		return new int[]{maxRank1, maxRank2};
	}
	
	public int[] getMaxSeries(List[] peaks, Spectrum s, int chargeDiff){
		Map<String, Peak> ionsPresented = new HashMap();
		LabelledPeak p;
		Set<Peak> matchedPeaks = new HashSet();
		if(peaks[0].size() == 0){
			return null;
		}
		String pep = ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide();
		int charge = ((LabelledPeak)peaks[0].get(0)).getPep().getCharge();
		for(int i = 0; i < peaks[0].size(); i++){
			p = (LabelledPeak)peaks[0].get(i);
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){
				ionsPresented.put(p.getType()+p.getPos()+"@"+p.getCharge(), (Peak)peaks[1].get(i));
				matchedPeaks.add((Peak)peaks[1].get(i));
			}
		}		
		
		int maxBSeries=0, maxYSeries=0, maxC=0;
		boolean prevB=false, prevY=false;
		for(int c = 1; c <= charge-chargeDiff; c++){
			int bseries=1, yseries=1;
			for(int i = 1; i <= pep.length(); i++){
				boolean currentB=false, currentY=false;
				for(int j = c; j <= c+chargeDiff; j++){
						currentB = currentB | ionsPresented.containsKey("b"+i+"@"+j);
				}
				
				if(prevB & currentB){
					bseries++;
				}
				if(!currentB){
					maxC = bseries > maxBSeries ? c : maxC;
					maxBSeries = bseries > maxBSeries ? bseries : maxBSeries;
					bseries=1;
				}
				prevB=currentB;
				
				for(int j = c; j <= c+chargeDiff; j++){
					currentY = currentY | ionsPresented.containsKey("y"+i+"@"+j);
				}
				if(prevY & currentY){ 
					yseries++;
				}
				if(!currentY){
					maxYSeries = yseries > maxYSeries ? yseries : maxYSeries;
					yseries=1;
				}
				prevY=currentY;
			}
			maxC = bseries > maxBSeries ? c : maxC;
			maxBSeries = bseries > maxBSeries ? bseries : maxBSeries;
			maxYSeries = yseries > maxYSeries ? yseries : maxYSeries;
		}
		//System.out.println("maxC is: " + maxC);
		return new int[]{maxBSeries, maxYSeries};
	}
	
	public void printIonsStatDetail(List[] peaks, Spectrum s){
		Map<String, Peak> ionsPresented = new HashMap();
		LabelledPeak p;
		if(peaks[0].size() == 0){
			return;
		}
		String pep = ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide();
		for(int i = 0; i < peaks[0].size(); i++){
			p = (LabelledPeak)peaks[0].get(i);
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){
				ionsPresented.put(p.getType()+p.getPos(), p);
			}
		}
		
		String ionType;
		for(int i = this.suffixIons.length-1; i >=0; i--){
			ionType = this.suffixIons[i];
			System.out.print(s.spectrumName + "\t" + this.padSpace(ionType, 9, false) + ": ");
			for(int j = pep.length(); j > 0; j--){
				if(ionsPresented.containsKey(ionType+j)){
					System.out.print("X");
//					p = (LabelledPeak)ionsPresented.get(ionType+j);
//					System.out.print(p.getCharge());
				}else{
					System.out.print("-");
				}
			}
			System.out.println();
		}
		System.out.println(s.spectrumName +  "\tPeptide  : " + pep);
		for(int i = 0; i < this.prefixIons.length; i++){
			ionType = this.prefixIons[i];
			System.out.print(s.spectrumName + "\t" + this.padSpace(ionType, 9, false) + ": ");
			for(int j = 0; j < pep.length(); j++){
				if(ionsPresented.containsKey(ionType+(j+1))){
					System.out.print("X");
//					p = (LabelledPeak)ionsPresented.get(ionType+(j+1));
//					System.out.print(p.getCharge());
				}else{
					System.out.print("-");
				}
			}
			System.out.println();
		}
		System.out.println();	
	}
	
	public void printIonsStatDetailWithCharge(List[] peaks, Spectrum s){
		Map<String, Peak> ionsPresented = new HashMap();
		LabelledPeak p;
		Set<Peak> matchedPeaks = new HashSet();
		if(peaks[0].size() == 0){
			return;
		}
		String pep = ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide();
		int charge = ((LabelledPeak)peaks[0].get(0)).getPep().getCharge();
		for(int i = 0; i < peaks[0].size(); i++){
			p = (LabelledPeak)peaks[0].get(i);
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){
				ionsPresented.put(p.getType()+p.getPos()+"@"+p.getCharge(), (Peak)peaks[1].get(i));
				matchedPeaks.add((Peak)peaks[1].get(i));
			}
		}
				
		String ionType;
		for(int i = this.suffixIons.length-1; i >=0; i--){
			ionType = this.suffixIons[i];
			System.out.print(s.spectrumName + "\t" + this.padSpace(ionType, 9, false) + ": ");
			List<Peak> presentedPeaks = new ArrayList<Peak>();
			for(int c = 1; c <= charge; c++){
				for(int j = pep.length(); j > 0; j--){
					if(ionsPresented.containsKey(ionType+j+"@"+c)){
						presentedPeaks.add(ionsPresented.get(ionType+(j)+"@"+c));
						System.out.print("X");
//											p = (LabelledPeak)ionsPresented.get(ionType+j);
//											System.out.print(p.getCharge());
					}else{
						System.out.print("-");
					}
				}
				System.out.print("     ");
			}
			for(int k = 0; k < presentedPeaks.size(); k++){
				Peak current = presentedPeaks.get(k);
				System.out.printf("[%.3f, %d]", current.getMass(), current.getRank());
			}
			System.out.println();
		}
		
		System.out.println(s.spectrumName +  "\tPeptide  : " + pep);
		for(int i = 0; i < this.prefixIons.length; i++){
			ionType = this.prefixIons[i];
			System.out.print(s.spectrumName + "\t" + this.padSpace(ionType, 9, false) + ": ");
			List<Peak> presentedPeaks = new ArrayList<Peak>();
			for(int c = 1; c <= charge; c++){
				for(int j = 0; j < pep.length(); j++){
					if(ionsPresented.containsKey(ionType+(j+1)+"@"+c)){
						System.out.print("X");
						presentedPeaks.add(ionsPresented.get(ionType+(j+1)+"@"+c));
//						p = (LabelledPeak)ionsPresented.get(ionType+j);
//						System.out.print(p.getCharge());
					}else{
						System.out.print("-");
					}
				}
				System.out.print("     ");
			}
			for(int k = 0; k < presentedPeaks.size(); k++){
				Peak current = presentedPeaks.get(k);
				System.out.printf("[%.3f, %d]", current.getMass(), current.getRank());
			}
			System.out.println();
		}
		System.out.println();
		List<Peak> toBePrint = new ArrayList();
		toBePrint.addAll(matchedPeaks);
		Collections.sort(toBePrint, PeakMassComparator.comparator);
		for(int i = 0; i < toBePrint.size(); i++){
			//System.out.println(toBePrint.get(i));
		}
		System.out.println(s.getPeptide() + " matched peaks: " + toBePrint.size());
	}
	
	public String padSpace(String orig, int length, boolean padFront){
		StringBuffer sb = new StringBuffer(orig);
		int diff = length - orig.length();
		if(padFront){
			for(int i = 0; i < diff; i++){
				sb.insert(0, " ");
			}	
		}else{
			for(int i = 0; i < diff; i++){
				sb.append(" ");
			}
		}
		return sb.toString();
	}
	
	//try to match peaks from two spectrum
	//return matched peaks in pair (theoretical, experimental) format
	public List[] matchSpectrum(Spectrum actual, double tolerance){
		List<Peak> list1 = new ArrayList<Peak>();
		List<Peak> list2 = new ArrayList<Peak>();
		Spectrum theoretical = this;
		List<Peak> pList1 = theoretical.getPeak(); 
		List<Peak> pList2 = actual.getPeak();
		int i = 0, j = 0, p = 0;
		double m1, m2, diff;
		while(i < pList1.size() && j < pList2.size()){
			m1 = pList1.get(i).getMass();
			m2 = pList2.get(j).getMass();
			
			if(Math.abs(m1 - m2) < tolerance){
				p = j;
				while(j < pList2.size() && Math.abs(m1 - pList2.get(j).getMass()) < tolerance){
					//System.out.println(actual.spectrumName + "\tmatched-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j).getIntensity());
					list1.add(pList1.get(i));
					list2.add(pList2.get(j));
					//once we find a match try match subsequent isotopic variance peaks
					int l = 1;
					for(int k = 1; k <= MAX_PROFILE_SPAN && j+k < pList2.size(); k++){
						m2 = pList2.get(j+k).getMass();
						m1 = m1 + 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
						diff = pList2.get(j+k).getMass() - pList2.get(j+k-1).getMass() - 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
						if(Math.abs(diff) < tolerance){
							//System.out.println(actual.spectrumName + "\tmatched-iso-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j+k).getIntensity());
							//list1.add(pList1.get(i));
							//list2.add(pList2.get(j+k));
							l++;
						}else{
							break;
						}
					}
					j+=l;
				}
				i++;
				j = p;
				continue;
			}
			if(m1 < m2){
				//System.out.println(actual.spectrumName + "\tunmatched-peaks:\t" + pList1.get(i) + "\t~\t");
				i++;

			}
			if(m2 < m1){
				//System.out.println(actual.spectrumName + "\tunmatched-peaks:\t" +  "\t\t\t\t\t\t\t~\t" + pList2.get(j).getMass() + "\t" + pList2.get(j).getIntensity());
				j++;
			}
			
		}
		return new List[]{list1, list2};
	}
	
	public SimpleGraph matchSpectrum2(Spectrum actual, double tolerance){
		List<Peak> list1 = new ArrayList<Peak>();
		List<Peak> list2 = new ArrayList<Peak>();
		Spectrum theoretical = this;
		List<Peak> pList1 = theoretical.getPeak(); 
		List<Peak> pList2 = actual.getPeak();
		SimpleGraph<Peak, DefaultEdge> matchingGraph = new SimpleGraph<Peak, DefaultEdge>(DefaultEdge.class);
		int i = 0, j = 0, p = 0;
		double m1, m2, diff;
		for(int c = 0; c < actual.getPeak().size();c++){
			matchingGraph.addVertex(actual.getPeak().get(c));
		}
		
		for(int c = 0; c < this.getPeak().size(); c++){
			matchingGraph.addVertex(this.getPeak().get(c));
		}
		
		while(i < pList1.size() && j < pList2.size()){
			m1 = pList1.get(i).getMass();
			m2 = pList2.get(j).getMass();
			
			if(Math.abs(m1 - m2) < tolerance){
				p = j;
				while(j < pList2.size() && Math.abs(m1 - pList2.get(j).getMass()) < tolerance){
					//System.out.println(actual.spectrumName + "\tmatched-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j).getIntensity());
					matchingGraph.addEdge(pList1.get(i), pList2.get(j));
					//once we find a match try match subsequent isotopic variance peaks
					int l = 1;
					for(int k = 1; k <= MAX_PROFILE_SPAN && j+k < pList2.size(); k++){
						m2 = pList2.get(j+k).getMass();
						m1 = m1 + 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
						diff = pList2.get(j+k).getMass() - pList2.get(j+k-1).getMass() - 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
						if(Math.abs(diff) < tolerance){
							//System.out.println(actual.spectrumName + "\tmatched-iso-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j+k).getIntensity());
							matchingGraph.addEdge(pList1.get(i), pList2.get(j+k));
							l++;
						}else{
							break;
						}
					}
					j+=l;
				}
				i++;
				j = p;
				continue;
			}
			if(m1 < m2){
				//System.out.println(actual.spectrumName + "\tunmatched-peaks:\t" + pList1.get(i) + "\t~\t");
				i++;

			}
			if(m2 < m1){
				//System.out.println(actual.spectrumName + "\tunmatched-peaks:\t" +  "\t\t\t\t\t\t\t~\t" + pList2.get(j).getMass() + "\t" + pList2.get(j).getIntensity());
				j++;
			}
		}
		
		return matchingGraph;
	}
	
	//fast matching with a slightly lower coverage of possibly matching peaks
	public SimpleGraph matchSpectrumx2(Spectrum actual, double tolerance){
		List<Peak> list1 = new ArrayList<Peak>();
		List<Peak> list2 = new ArrayList<Peak>();
		Spectrum theoretical = this;
		List<Peak> pList1 = theoretical.getPeak(); 
		List<Peak> pList2 = actual.getPeak();
		SimpleGraph<Peak, DefaultEdge> matchingGraph = new SimpleGraph<Peak, DefaultEdge>(DefaultEdge.class);
		int i = 0, j = 0, p = 0;
		double m1, m2, diff;
		Peak p1, p2;
		for(int c = 0; c < actual.getPeak().size();c++){
			matchingGraph.addVertex(actual.getPeak().get(c));
		}
		
		pList2.add(new Peak(100000, 0)); //we add a end-guard peak to end of the list
		                                 //makes subsquent iteration easier
		for(int c = 0; c < this.getPeak().size(); c++){
			matchingGraph.addVertex(this.getPeak().get(c));
		}
		p1 = pList1.get(0);
		p2 = pList2.get(0);
		for(i = 0; i < pList1.size(); i++){
			p1 = pList1.get(i);
			for(; j < pList2.size(); j++){
				p2 = pList2.get(j);
				//scanning left to reach lower boundary for a theo-peak
				if(p1.getMass() - p2.getMass() > tolerance){
					continue;
				//once we surpass lower boundary, we check whether we are over 
				//upper boundary if it is we skip this theo-peak
				}else if(p2.getMass() - p1.getMass() > tolerance){
					break;
				//now we are within tolerance boundary	
				}else{
					matchingGraph.addEdge(p1, p2);
//					p = j;
//					do{
//						matchingGraph.addEdge(p1, p2);
//						p2 = pList2.get(++p);
//					}while(p2.getMass() - p1.getMass() < tolerance);
//					for(int k = 1; k <= MAX_PROFILE_SPAN && j+k < pList2.size(); k++){
//						diff = pList2.get(j+k).getMass() - pList2.get(j+k-1).getMass() - 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
//						if(Math.abs(diff) < tolerance){
//							//System.out.println(actual.spectrumName + "\tmatched-iso-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j+k).getIntensity());
//							matchingGraph.addEdge(pList1.get(i), pList2.get(j+k));
//						}else{
//							break;
//						}
//					}
				}
				
			}
		}
		pList2.remove(pList2.size()-1);
		//System.out.println("there are total of " + matchingGraph.edgeSet().size() + " mathcing within tolerance");
		return matchingGraph;
	}
	public SimpleMatchingGraph getMatchGraph(Spectrum actual, double tolerance){
		return getMatchGraph(actual, tolerance, false);
	}
	
	public SimpleMatchingGraph getMatchGraph(Spectrum actual, double tolerance, boolean detail){
		return getMatchGraph(actual, tolerance, detail, 0);
	}
	
	//Note: to speed things up, if there is less than N matched peaks in the spectrum we 
	//return an empty match graph instead
	public SimpleMatchingGraph getMatchGraph(Spectrum actual, double tolerance, boolean detail, int minMatchPeaks){
//		List<Peak> list1 = new ArrayList<Peak>(200);
//		List<Peak> list2 = new ArrayList<Peak>(200);
		Spectrum theoretical = this;
		List<Peak> pList1 = theoretical.getPeak(); 
		List<Peak> pList2 = actual.getPeak();
		SimpleMatchingGraph matchingGraph = new SimpleMatchingGraph();
		pList2.add(new Peak(100000, 0)); //we add a end-guard peak to end of the list
        //makes subsquent iteration easier
		//SimpleMatchingGraph matchingGraph = SimpleMatchGraphFactory.createSimpleMatchGraph();
		int i = 0, j = 0, p = 0;
		double m1, m2, diff;
		Peak p1, p2;
		p1 = pList1.get(0);
		p2 = pList2.get(0);
		for(i = 0; i < pList1.size(); i++){
			p1 = pList1.get(i);
			if(detail)
				System.out.println("matching " + p1 + " starting " + j);
			for(; j < pList2.size(); j++){
				p2 = pList2.get(j);
				//scanning left to reach lower boundary for a theo-peak
				if(p1.getMass() - p2.getMass() > tolerance){
					if(detail){
						System.out.println("not-matching: " + p2 + "~");
					}
					continue;
				//once we surpass lower boundary, we check whether we are over 
				//upper boundary if it is we skip this theo-peak
				}else if(p2.getMass() - p1.getMass() > tolerance){
					if(detail){
						System.out.println("not-matching: " + p2 + "~");
					}
					break;
				//now we are within tolerance boundary	
				}else{
					//matchingGraph.addEdge(p2, p1);
					p = j;
					do{
						matchingGraph.addVertex(p2, 1);
						matchingGraph.addVertex(p1, 2);
						matchingGraph.addEdge(p2, p1);
						//detail = true;
						if(detail){
							System.out.println("matching: " + p2 + " ~ " + p1);
						}
						//detail = false;
						p2 = pList2.get(++p);
					}while(p2.getMass() - p1.getMass() < tolerance);
					break;
//					for(int k = 1; k <= MAX_PROFILE_SPAN && j+k < pList2.size(); k++){
//						diff = pList2.get(j+k).getMass() - pList2.get(j+k-1).getMass() - 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
//						if(Math.abs(diff) < tolerance){
//							//System.out.println(actual.spectrumName + "\tmatched-iso-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j+k).getIntensity());
//							matchingGraph.addEdge(pList2.get(j+k), pList1.get(i));
//						}else{
//							break;
//						}
//					}
				}	
			}
		}
		pList2.remove(pList2.size()-1);
		//System.out.println("matching graph size: " + matchingGraph.vertexSet(1).size() + "\t" + matchingGraph.vertexSet(2).size());
		if(matchingGraph.vertexSet(SimpleMatchingGraph.Observed).size() < minMatchPeaks){
			return SimpleMatchGraphFactory.getEmptyGraph();
		}
		//System.out.println("matching graph size: " + matchingGraph.vertexSet(SimpleMatchingGraph.Observed).size());

		for(int c = 0; c < pList2.size();c++){
			matchingGraph.addVertex(actual.getPeak().get(c), 1);
		}
		
		
		for(int c = 0; c < pList1.size(); c++){
			matchingGraph.addVertex(pList1.get(c), 2);
		}

		//System.out.println("there are total of " + matchingGraph.edgeSet().size() + " mathcing within tolerance");
		return matchingGraph;

	}
	
	public SimpleMatchingGraph getPPMMatchGraph(Spectrum actual, double ppm, boolean detail, int minMatchPeaks){
//		List<Peak> list1 = new ArrayList<Peak>(200);
//		List<Peak> list2 = new ArrayList<Peak>(200);
		Spectrum theoretical = this;
		List<Peak> pList1 = theoretical.getPeak(); 
		List<Peak> pList2 = actual.getPeak();
		SimpleMatchingGraph matchingGraph = new SimpleMatchingGraph();
		pList2.add(new Peak(100000, 0)); //we add a end-guard peak to end of the list
        //makes subsquent iteration easier
		//SimpleMatchingGraph matchingGraph = SimpleMatchGraphFactory.createSimpleMatchGraph();
		int i = 0, j = 0, p = 0;
		double m1, m2, diff;
		Peak p1, p2;
		p1 = pList1.get(0);
		p2 = pList2.get(0);
		for(i = 0; i < pList1.size(); i++){
			p1 = pList1.get(i);
			if(detail)
				System.out.println("matching " + p1 + " starting " + j);
			for(; j < pList2.size(); j++){
				p2 = pList2.get(j);
				//scanning left to reach lower boundary for a theo-peak
				//System.out.println("error: " + (p1.getMass() - p2.getMass())*1000000/p2.getMass() + "\t" + p1.getMass() + "\t" + p2.getMass());
				//System.out.println("error: " + (p2.getMass() - p1.getMass())*1000000/p2.getMass() + "\t" + p1.getMass() + "\t" + p2.getMass());
				if((p1.getMass() - p2.getMass())*1000000/p2.getMass() > ppm){
					if(detail){
						System.out.println("not-matching: " + p2 + "~");
					}
					continue;
				//once we surpass lower boundary, we check whether we are over 
				//upper boundary if it is we skip this theo-peak
				}else if((p2.getMass() - p1.getMass())*1000000/p2.getMass() > ppm){
					if(detail){
						System.out.println("not-matching: " + p2 + "~");
					}
					break;
				//now we are within tolerance boundary	
				}else{
					//matchingGraph.addEdge(p2, p1);
					p = j;
					do{
						matchingGraph.addVertex(p2, 1);
						matchingGraph.addVertex(p1, 2);
						matchingGraph.addEdge(p2, p1);
						if(detail){
							System.out.println("matching: " + p2 + " ~ " + p1);
						}
						p2 = pList2.get(++p);
					}while((p2.getMass() - p1.getMass())*1000000/p2.getMass() < ppm);
					break;
//					for(int k = 1; k <= MAX_PROFILE_SPAN && j+k < pList2.size(); k++){
//						diff = pList2.get(j+k).getMass() - pList2.get(j+k-1).getMass() - 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
//						if(Math.abs(diff) < tolerance){
//							//System.out.println(actual.spectrumName + "\tmatched-iso-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j+k).getIntensity());
//							matchingGraph.addEdge(pList2.get(j+k), pList1.get(i));
//						}else{
//							break;
//						}
//					}
				}	
			}
		}
		pList2.remove(pList2.size()-1);
		//System.out.println("matching graph size: " + matchingGraph.vertexSet(1).size() + "\t" + matchingGraph.vertexSet(2).size());
		if(matchingGraph.vertexSet(SimpleMatchingGraph.Observed).size() < minMatchPeaks){
			return SimpleMatchGraphFactory.getEmptyGraph();
		}
		//System.out.println("matching graph size: " + matchingGraph.vertexSet(SimpleMatchingGraph.Observed).size());

		for(int c = 0; c < pList2.size();c++){
			matchingGraph.addVertex(actual.getPeak().get(c), 1);
		}
		
		
		for(int c = 0; c < pList1.size(); c++){
			matchingGraph.addVertex(pList1.get(c), 2);
		}

		//System.out.println("there are total of " + matchingGraph.edgeSet().size() + " mathcing within tolerance");
		return matchingGraph;

	}
	
	public SimpleMatchingGraph getMatchGraph(ArraySpectrum actual, double tolerance, boolean detail, int minMatchPeaks){
		ArraySpectrum theoretical = (ArraySpectrum)this;
		Peak[] pList1 = theoretical.intensities;
		Peak[] pList2 = new Peak[actual.intensities.length+1];
		double[] masses1 = theoretical.masses;
		double[] masses2 = new double[actual.masses.length+1];
		
		for(int i = 0; i < actual.masses.length; i++){
			masses2[i] = actual.masses[i];
			pList2[i] = actual.intensities[i];
		}
		
		SimpleMatchingGraph matchingGraph = new SimpleMatchingGraph();
		pList2[pList2.length-1]= new Peak(100000, 0); //we add a end-guard peak to end of the list
        masses2[masses2.length-1]=1000000;
		//makes subsquent iteration easier
		//SimpleMatchingGraph matchingGraph = SimpleMatchGraphFactory.createSimpleMatchGraph();
		int i = 0, j = 0, p = 0;
		double m1, m2, diff;
		Peak p1, p2;
		p1 = pList1[0];
		p2 = pList2[0];
		int length1 = pList1.length;
		int length2 = pList2.length;
		for(i = 0; i < length1; i++){
			p1 = pList1[i];
			if(detail)
				System.out.println("matching " + p1 + " starting " + j);
			for(; j < length2; j++){
				p2 = pList2[j];
				//scanning left to reach lower boundary for a theo-peak
				if(masses1[i] - masses2[j] > tolerance){
					continue;
				//once we surpass lower boundary, we check whether we are over 
				//upper boundary if it is we skip this theo-peak
				}else if(masses2[j] - masses1[i] > tolerance){
					break;
				//now we are within tolerance boundary	
				}else{
					//matchingGraph.addEdge(p2, p1);
					p = j;
					do{
						matchingGraph.addVertex(p2, 1);
						matchingGraph.addVertex(p1, 2);
						matchingGraph.addEdge(p2, p1);
						if(detail){
							System.out.println("matching: " + p2 + " ~ " + p1);
						}
						p2 = pList2[++p];
					}while(masses2[p] - masses1[i] < tolerance);
					break;
//					for(int k = 1; k <= MAX_PROFILE_SPAN && j+k < pList2.size(); k++){
//						diff = pList2.get(j+k).getMass() - pList2.get(j+k-1).getMass() - 1.0/(double)((LabelledPeak)pList1.get(i)).getCharge();
//						if(Math.abs(diff) < tolerance){
//							//System.out.println(actual.spectrumName + "\tmatched-iso-peaks\t" + pList1.get(i) + "\t~\t" + m2 + "\t" + pList2.get(j+k).getIntensity());
//							matchingGraph.addEdge(pList2.get(j+k), pList1.get(i));
//						}else{
//							break;
//						}
//					}
				}	
			}
		}
		
		//System.out.println("matching graph size: " + matchingGraph.vertexSet(1).size() + "\t" + matchingGraph.vertexSet(2).size());
		if(matchingGraph.vertexSet(SimpleMatchingGraph.Observed).size() < minMatchPeaks){
			return SimpleMatchGraphFactory.getEmptyGraph();
		}
		//System.out.println("matching graph size: " + matchingGraph.vertexSet(SimpleMatchingGraph.Observed).size());
		for(int c = 0; c < length2;c++){
			matchingGraph.addVertex(pList2[c], 1);
		}
		
		
		for(int c = 0; c < length1; c++){
			matchingGraph.addVertex(pList1[c], 2);
		}

		//System.out.println("there are total of " + matchingGraph.edgeSet().size() + " mathcing within tolerance");
		return matchingGraph;
	}
	
	
	
	public List[] refineMatchedSpectrum(List[] matchedList, Spectrum s){
		Map<String, Peak> primaryPresented = new HashMap();
		LabelledPeak p;
		List refined1 = new ArrayList(), refined2 = new ArrayList();
		for(int i = matchedList[0].size()-1; i >= 0; i--){
			p = (LabelledPeak)matchedList[0].get(i);
			if(p.isPrimary()){
//				System.out.println("adding primary ions: " + p.getType() + p.getPos() + "@" + p.getCharge());
				primaryPresented.put(p.getPep().getPeptide()+p.getType() + p.getPos() + "@" + p.getCharge(), p);
				refined1.add(p);                      //adding primary ions to the list
				refined2.add(matchedList[1].get(i));
			}else{
				if(primaryPresented.containsKey(p.getPep().getPeptide()+p.getType().substring(0,1) + p.getPos() + "@" + p.getCharge())){
//					System.out.println("seeing secondary ions: " + p.getType().substring(0,1) + p.getPos() + "@" + p.getCharge());
					refined1.add(p);                     //we only added secondary ions (e.g -H2O, -NH3 etc..) only
					refined2.add(matchedList[1].get(i)); //if the primary ions are presented
				}
			}
		}
		return new List[]{refined1, refined2};
	}
	
	public SimpleMatchingGraph refineMatchedSpectrum(SimpleMatchingGraph mGraph, Spectrum s){
		Map<String, Peak> primaryPresented = new HashMap();
		LabelledPeak p;
		Iterator<LabelledPeak> it = mGraph.vertexSet(2).iterator();
		while(it.hasNext()){
			p = it.next();
			if(p.isPrimary() && mGraph.getNeighbors(p).size() > 0){
				primaryPresented.put(p.getPep().getPeptide()+p.getType() + p.getPos() + "@" + p.getCharge(), p);
				//System.out.println("adding primary ions: " + p.getPep().getPeptide()+p.getType() + p.getPos() + "@" + p.getCharge());
			}
		}
		
		it = mGraph.vertexSet(2).iterator();
		while(it.hasNext()){
			p = it.next();
			//System.out.println("checking ions: " + p.getPep().getPeptide()+p.getType().substring(0,1) + p.getPos() + "@" + p.getCharge());
			if(!primaryPresented.containsKey(p.getPep().getPeptide()+p.getType().substring(0,1) + p.getPos() + "@" + p.getCharge())){
				//System.out.println("removing only 2nd ions match");
				mGraph.removeAllNeighbor(p);
			}
		}
		return mGraph;
	}
	
	private void splitPeaks(List[] peaks, List[] l1, List[] l2){
		String firstPeptide = ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide();
		splitPeaks(peaks, l1, l2, firstPeptide, null);
	}
	
	
	//note 2nd peptide can be null argument, be caution when use that
	private void splitPeaks(List[] peaks, List[] l1, List[] l2, String peptide1, String peptide2){
		if(peaks[0].size() == 0){
			return;
		}
		String firstPeptide = peptide1;
		LabelledPeak p;
		Peak p2;
		for(int i = 0; i < peaks[0].size(); i++){
			p = (LabelledPeak)peaks[0].get(i);
			p2 = (Peak)peaks[1].get(i);
			//System.out.println("matching " +p.getPep().getPeptide() + " to " + firstPeptide);
			if(p.getPep().getPeptide().equals(firstPeptide)){
				l1[0].add(p);
				l1[1].add(p2);
			}else{
				l2[0].add(p);
				l2[1].add(p2);

			}
		}
		//System.out.println("split peaks count: " + peaks[0].size() + "\t" + peaks[1].size());
	}
	
	public static boolean isLinkedPeak(Peptide pep, LabelledPeak p){
		boolean prefixLink = p.isPrefixPeak() && p.getPos() >= pep.getLinkedPos();
		boolean suffixLink = p.isSuffixPeak() && p.getPos() > pep.getPeptide().length() - pep.getLinkedPos();
		return prefixLink || suffixLink;
		//return false;
	}
	
	//shift all peaks at certain position by certain amount
	//peak list must be sorted by the mass-to-charge ratio
	public void shiftPeaks(List<LabelledPeak> pList, int[] pos, double[] mass, int length){
		System.out.println("shifting at position: "  + pos[0]);
		LabelledPeak current;
		int l = 0;
//		int length = pList.get(pList.size()-1).getPos(); //lenght of the peptide
		double mod = 0.0;
		for(int i = 0; i < pList.size(); i++){
			current = pList.get(i);
			for(int j = 0; j < pos.length; j++){
				if(( (current.getType().contains("b") || current.getType().contains("a"))
						&& current.getPos() >= pos[j] )
						|| (current.getType().contains("y") && current.getPos() > length-pos[j]) ){
					current.setMoz(current.getMass()+mass[j]/(double)current.getCharge());
				}else{
					break;
				}
			}
			
		}
		Collections.sort(pList, new PeakMassComparator());
	}
	
	public void addChargedPeaks(List<LabelledPeak> pList, int[] pos, int charge){
		//System.out.println("before adding we have: " + pList.size());
		//System.out.println("adding charge up to " + charge);
		LabelledPeak current, charged;
		int l = 0;
		int length = pList.get(pList.size()-1).getPos(); //lenght of the peptide
		//System.out.println("peptide length is: " + length);
		double mod = 0.0;
		List<LabelledPeak> added = new Vector();
		int count = 0;
		List toBeRemoved = new ArrayList();
		Set toBeRemoved2 = new HashSet();
		for(int i = 0; i < pList.size(); i++){
			current = pList.get(i);
			//add charged peaks for all cross-linked species
			for(int j = 0; j < pos.length; j++){
				if(((current.getType().contains("b") || current.getType().contains("a"))&& current.getPos() >= pos[j])
						|| (current.getType().contains("y") && current.getPos() > length-pos[j])){
					if(current.getCharge() == 1){
						count++;
						//for(int c = current.getCharge()+1; c <= charge; c++){
						//int minCharge=charge-2;
						int minCharge = 1;
						minCharge = minCharge <= 0 ? 1 : minCharge;
						int maxCharge = charge;
						if(charge > 5){
							//maxCharge = 5;
						}
						maxCharge = charge;
						for(int c = minCharge; c <= maxCharge; c++){
							charged = new LabelledPeak((current.getMass() + Mass.PROTON_MASS*(c-1))/(double)c,
									current.getIntensity(), 
									current.getType(),
									current.getPos(),
									(short)c);
							added.add(charged);
						}
					}
					toBeRemoved2.add(current);
				}else{
					break;
				}
			}
		}
		pList.addAll(added);
		pList.removeAll(toBeRemoved2);
		Collections.sort(pList, new PeakMassComparator());
		//System.out.println("adding charged-1 ions: " + count);
		//System.out.println("after adding we have: " + pList.size());
		LabelledPeak p1, p2;
		for(int i = 1; i < pList.size(); i++){
			p1 = pList.get(i-1);
			p2 = pList.get(i);
			if(Math.abs(p1.getMass() - p2.getMass()) < 0.0000000001 
					&& p1.getType() == p2.getType()
					&& p1.getCharge() == p2.getCharge()
					&& !p1.equals(p2)){
				//System.out.println(p1 +  " and " + p2 + " is same peak");
				toBeRemoved.add(p1);
			}
		}
		pList.removeAll(toBeRemoved);
		//System.out.println("after removal we have " + pList.size());
	}
	
	public double[] analyzeAnnotation(Spectrum s, String p){
		return analyzeAnnotation(s, p, 0.3);
	}
	
	public double[] analyzeAnnotation(Spectrum s, String p, double tolerance){
		return analyzeAnnotation(s, p, tolerance, true);
	}
	
	public double[] analyzeAnnotation(Spectrum s, String p, double tolerance, boolean detail){
		List[] matching = this.matchSpectrum(s, tolerance);
		double totalFraction = this.explainedPeaks2(matching, s, tolerance, detail);
		//this.splitPeaks(matching, l1, l2);
		double[] b_y_fraction = this.IonsStat(matching, s, detail);
		if(detail)	this.printIonsStatDetailWithCharge(matching, s);
		double merror = avgError(matching);
		double intensity = intensity(matching);
		//this.printIonsStatDetail(l1, s);
		//System.out.println();
		//this.printIonsStatDetail(l2, s);
		double[] stat = new double[6];
		stat[0] = totalFraction;
		stat[1] = b_y_fraction[0];
		stat[2] = b_y_fraction[1];
		stat[3] = b_y_fraction[2];
		stat[4] = b_y_fraction[3];
		stat[5] = merror;
		return stat;
	}
		
	public double[] analyzeMixtureAnnotation(Spectrum s, String p1, String p2){
		return analyzeMixtureAnnotation(s, p1, p2, TheoreticalSpectrum.MS2Tolerance);
	}
	public double[] analyzeMixtureAnnotation(Spectrum s, String p1, String p2, double tolerance){
		return analyzeMixtureAnnotation(s, p1, p2, tolerance, true);
	}
	public double[] analyzeMixtureAnnotation(Spectrum s, String p1, String p2, double tolerance, boolean detail){
		String[] peptides = this.peptide.split(" & ");
		peptides[0] = peptides[0].replaceAll("[0-9\\.\\+\\-]", "");
		peptides[1] = peptides[1].replaceAll("[0-9\\.\\+\\-]", "");
		if(detail)
			System.out.println("peptides are: " + peptides[0] + " and " + peptides[1]);		
		if(peptides.length != 2){
			System.err.println("warnining: mixture spectrum peptides is not instatiate correctly");
		}
		List[] matching = this.matchSpectrum(s, tolerance);
		List[] l1 = new List[] {new ArrayList(), new ArrayList()};
		List[] l2 = new List[] {new ArrayList(), new ArrayList()};
		//make sure the order is the same
		
		double totalFraction = this.explainedPeaks2(matching, s, tolerance, detail);
		//System.out.println("peptides " + p1);
		this.splitPeaks(matching, l1, l2, peptides[0].split("\\.")[0], peptides[1].split("\\.")[0]);
		//this.splitPeaks(matching, l1, l2);
//		System.out.println("combine size: " + matching[0].size() + "\t" + matching[1].size());
//		System.out.println("1size: " + l1[0].size()+ "\t" + l1[1].size());
//		System.out.println("2size: " + l2[0].size() + "\t" + l2[1].size());
		double[] b_y_fraction1 = this.IonsStat(l1, s, detail);
		double[] b_y_fraction2 = this.IonsStat(l2, s, detail);
		double total = s.sumOfPeaks()*s.sumOfPeaks();
		double merror1 = avgError(l1);
		double merror2 = avgError(l2);
		//this.printIonsStatDetail(l1, s);
		//System.out.println();
		//this.printIonsStatDetail(l2, s);
		if(detail){
			this.printIonsStatDetailWithCharge(l1, s);
			System.out.println();
			this.printIonsStatDetailWithCharge(l2, s);
		}
		double[] stat = new double[17];
		stat[0] = totalFraction;
		stat[1] = b_y_fraction1[0];
		stat[2] = b_y_fraction1[1];
		stat[3] = b_y_fraction2[0];
		stat[4] = b_y_fraction2[1];
		stat[5] = b_y_fraction1[2];
		stat[6] = b_y_fraction1[3];
		stat[7] = b_y_fraction2[2];
		stat[8] = b_y_fraction2[3];
		stat[9] = merror1;
		stat[10] = merror2;
		stat[11] = this.explainedPeaks2(l1, s, tolerance, detail);
		stat[12] = this.explainedPeaks2(l2, s, tolerance, detail);
		stat[13] = b_y_fraction1[4];
		stat[14] = b_y_fraction1[5];
		stat[15] = b_y_fraction2[4];
		stat[16] = b_y_fraction2[5];
		return stat;
	}

	public double avgError(List[] peaks){
		LabelledPeak lp;
		Peak p;
		if(peaks[0].size()==0){
			return 0.0;
		}
		double massDiff = 0.0;
		for(int i = 0; i < peaks[0].size(); i++){
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){	
				lp = (LabelledPeak)peaks[0].get(i);
				p = (Peak)peaks[1].get(i);
				massDiff += Math.abs(lp.getMass() - p.getMass());
			}
		}
		return massDiff/peaks[0].size();
	}
	
	
	public double intensity(List[] peaks){
		double intensity = 0.0;
		for(int i = 0; i < peaks[0].size(); i++){
				Peak p = (Peak)peaks[1].get(i);
				intensity += p.getIntensity();
		}
		return intensity;
	}
	
	public SpectrumLib createTheoSpectrumLibrary(String peptidesFile){
		HashMap<String, List<Spectrum>> table = new HashMap();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(peptidesFile));
			String currentLine = bf.readLine(); //skip header
			String[] tokens;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			Vector v;
			while(currentLine != null){
				th = new TheoreticalSpectrum(currentLine);
				currentLine = bf.readLine();
				v = new Vector();
				v.add(th);
				table.put(th.peptide, v);
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
		return new SpectrumLib(table);
	}
	
	public void printSpectPlotAnnotationFile(Spectrum s, String file, double masstolerance){
		SimpleMatchingGraph g = this.getMatchGraph(s, masstolerance);
		printSpectPlotAnnotationFile(g, s, file, masstolerance);
	}
	
	public void printSpectPlotAnnotationFile(SimpleMatchingGraph g, Spectrum s, String file, double masstolerance){
		try{
			BufferedWriter bo = new BufferedWriter(new FileWriter(file));
			bo.append("BEGIN IONS\n");
			Map<Peak, String> colorMap = new HashMap();
			//int counter = 0;
			for(Iterator<Peak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
				Peak current = it.next();
				List<Peak> neighbors = g.getNeighbors(current);
				String[] peps=this.peptide.split(" & ");
				String color = null;
				String type = null;
				//System.out.println("current peaks: " + current);
				for(Iterator it2 = neighbors.iterator(); it2.hasNext();){
					//counter++;
					LabelledPeak current2 = (LabelledPeak)it2.next();
					//System.out.println("peptide is: " + current2.getPep() + "\t" + current2.getPep().toString().length());
					//System.out.println("current label: " + current2);
					if(current2.getPep().toString().equals(peps[0]) && current2.isPrimary()){
						if(color == null || color.equals("blue")){
							color = "blue";
						}else{
							color = "green";
							//break;
						}
						if(type != null){
							type = "";
						}else{
							type = type+ ",";
						}
						type=type+"\\alpha ";
						type = current2.getType();
						type = type.substring(0,1)+current2.getPos();
						if(current2.getCharge() > 1){
							type = type.substring(0,1)+"@"+"^{"+current2.getPos()+"}";
							type = type + "_{" + current2.getCharge()+"}";
						}
					}else if(current2.getPep().toString().equals(peps[1]) && current2.isPrimary()){
						if(color == null || color.equals("red")){
							color = "red";
						}else{
							color = "green";
							//break;
						}
						if(type != null){
							type = "";
						}else{
							type = type+ ",";
						}
						type=type+"\\beta ";
						type = current2.getType();
						type = type.substring(0,1)+current2.getPos();
						if(current2.getCharge() > 1){
							type = type.substring(0,1)+"@"+"^{"+current2.getPos()+"}";
							type = type + "_{" + current2.getCharge()+"}";
						}

					}		
				}
				if(neighbors.size() > 0 && color == null){
					color = "black";
					type = "";
				}
				if(neighbors.size() > 0){
					//System.out.println("printing peaks: " + current);
					colorMap.put(current, color);
					bo.append(current.getMass() + "\t" 
							+ current.getIntensity()+"\t"
							+ type + "\t\t"
							+ "1\t" + 1 +"\t"
							+ color +"\n");
				}
				//System.out.println("color is: " + color);
			}
			//try matching unfragmented precursors
			for(Iterator<Peak> it = s.getPeaks().iterator(); it.hasNext();){
				Peak current = it.next();
				if(Math.abs(current.getMass() - s.parentMass) < masstolerance){
					bo.append(current.getMass() + "\t" 
							+ current.getIntensity()+"\t"
							+ "Pr" + "\t\t"
							+ "1\t" + 1 +"\t"
							+ "black" +"\n");
				}
				if(s.charge > 0){
					if(Math.abs(current.getMass() - 
							(s.parentMass - Mass.WATER/s.charge)) < masstolerance
							|| Math.abs((s.parentMass - Mass.NH3/s.charge)) < masstolerance
							|| Math.abs(current.getMass() - 
									(s.parentMass - 2*Mass.WATER/s.charge)) < masstolerance
							|| Math.abs(current.getMass() - 
											(s.parentMass - (Mass.WATER+Mass.NH3)/s.charge)) < masstolerance){
						bo.append(current.getMass() + "\t" 
								+ current.getIntensity()+"\t"
								+ "Pr" + "\t\t"
								+ "1\t" + 1 +"\t"
								+ "black" +"\n");
					}
				}
				
			}
			//System.out.println("total matched peaks " + counter);
			for(Iterator<Peak> it = g.vertexSet(g.Theoretical).iterator(); it.hasNext();){
				LabelledPeak current2 = (LabelledPeak)it.next();
				List<Peak> neighbors = g.getNeighbors(current2);
				if(neighbors.size() > 0){
					double max=-1000;
					Peak highest =null;
					for(int i= 0; i < neighbors.size(); i++){
						Peak neigh = neighbors.get(i);
						highest = neigh.getIntensity() > max ? neigh : highest;
						max = neigh.getIntensity() > max ? neigh.getIntensity() : max;

					}
					
					for(int i = 0; i < neighbors.size(); i++){
						Peak neigh = neighbors.get(i);
						String type = current2.getType();
						type = type.substring(0,1)+current2.getPos();
						if(current2.getCharge() > 1){
							type = type.substring(0,1)+"@"+"^{"+current2.getPos()+"}";
							type = type + "_{" + current2.getCharge()+"}";
						}
						if(neigh != highest || current2.getType().length()!=1){
							type = "";
						}
						String color = colorMap.get(neigh);
//						bo.append(neigh.getMass() + "\t" 
//								+ neigh.getIntensity()+"\t"
//								+ type + "\t\t"
//								+ "1\t" + current2.getCharge()+"\t"
//								+ color +"\n");
					}
				}
			}
			bo.append("END IONS");
			bo.flush();
			bo.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		}
	}
	
	public static void testTheoreticalMasses(){
		String libFile = ".\\MSPLib\\Lib\\yeast.msp";
		String spectrumFile = "..\\mixture_linked\\yeast_specLib.mgf";
		String annotationFile ="..\\mixture_linked\\trps\\result.txt";
		TheoreticalSpectrum t;
		String[] prefixIons = {"b"};
		String[] suffixIons = {"y"};
		SpectrumLib lib1 = new SpectrumLib(libFile, "MSP");
		lib1.removeModSpectra();
		System.out.println("finish loading");
//		lib1.annoateSpectrumFromInspectFile(annotationFile);
		lib1.computeRank();
		lib1.toRelativeIntensity();
		List<Spectrum> subset = SpectrumUtil.getSpectra("..\\mixture_linked\\t", lib1);
		Map<String, String> map = FileIOUtils.createTableFromFile("..\\mixture_linked\\t2", 0, 0);
		Iterator<Spectrum> it = subset.listIterator();
		Iterator<String> it2 = map.keySet().iterator();
		RankBaseScoreLearner peakscorer = new RankBaseScoreLearner(lib1);
		peakscorer.getIonsCount();
		lib1 = null;
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		Spectrum s;
		List[] matchedpeaks;
		while(it.hasNext()){
			s = it.next();
			if(!s.peptide.contains("Scan") && /*s.peptide.contains("+") &&*/ (s.charge >= 1)){
				//s.filterPeaks(6*s.peptide.length());
				//s.windowFilterPeaks(6, 25); //inspect parameter
				t = new TheoreticalSpectrum(it2.next());
				matchedpeaks = t.matchSpectrum(s, 0.5);
				t.explainedPeaks2(matchedpeaks, s, 0.5);
				t.IonsStat(matchedpeaks, s);
				t.printIonsStatDetailWithCharge(matchedpeaks, s);
				double score = scorer.compare(t, s);
				System.out.println("score is: " + score);
				//return;
			}
		}
	}
	
	public static void testTheoreticalMixture(){
		String spectrumFile = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
		String mixtureFile = ".\\mixture_linked\\tk090204_WKarzai_Tryp.mgf";
		String annotationFile =".\\mixture_linked\\trps\\result.txt";
		String tripletFile =".\\mixture_linked\\triplet.txt";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
		System.out.println("Done loading spectrum library");
		lib1.annoateSpectrumFromInspectFile(annotationFile);
		SpectrumLib mix = new SpectrumLib(mixtureFile, "MGF");
		System.out.println("Done loading linked-mixture library");
		lib1.toRelativeIntensity();
		mix.toRelativeIntensity();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); //skip header
			String[] tokens;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				System.out.println("line is: " + currentLine);
				tokens = currentLine.split("\t");
				m = mix.getSpectra(tokens[0]).get(0);
				s1 = lib1.getSpectra(tokens[1]).get(0);
				s2 = lib1.getSpectra(tokens[2]).get(0);
				currentLine = bf.readLine();
				System.out.println("peptides are: " + s1.peptide + " & " + s2.peptide);
				//mix.filterPeaks(99);
				//th = new TheoreticalSpectrum(s1.peptide, s2.peptide);
				//th = new TheoreticalSpectrum(s1, s2);
				testMultipleLysPosition(s1, s2, m);
//				System.out.println("done creating theo-spect");
//				System.out.println(m.spectrumName + "(" + th.peptide + ") has explained b/y ions fraction: " + 
//						th.explainedPeaks2(th, m, 0.50));
				//return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void testTheoreticalMixture2(){
		String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		String mixSpectrumFile = "..\\mixture_linked\\yeast_data\\klc_122007p_yeast_digest1.mgf";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		lib1.removeModSpectra();
		lib1.computeRank();
		System.out.println("Done loading spectrum library");
		//SpectrumLib mix = lib1.createMix("..\\mixture_linked\\ecoli_mixture2.name", 400, 1, 1, 0, 1, 2000, false);
		SpectrumLib mix = new SpectrumLib(mixSpectrumFile, "MGF");
		System.out.println("Done loading linked-mixture library");
		String tripletFile = "..\\mixture_linked\\mixtureAnnotation.txt";
		tripletFile = "..\\mixture_linked\\temp";
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel2.txt";
		//lib1.toRelativeIntensity();
		mix.toRelativeIntensity();
//		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
//		SpectrumIonRankLearner learner = new SpectrumIonRankLearner(lib1);
//		PeakComparator peakscorer = learner.createComparatorSet();
		RankBaseScoreLearner peakscorer = new RankBaseScoreLearner(lib1);
		peakscorer.getIonsCount();
		//SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		RankBaseMixtureSpectrumScorer scorer = new RankBaseMixtureSpectrumScorer(peakscorer);
		SimpleMatchingGraph g;
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); 
			String[] tokens, tokens2;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				tokens = currentLine.split("\t");
				//tokens2 = tokens[0].split(" & ");
				m = mix.getSpectra(tokens[0]).get(0);
				m.windowFilterPeaks(15, 25);
				m.computePeakRank();				
				th = new TheoreticalSpectrum(tokens[1], tokens[2]);
				th.analyzeMixtureAnnotation(m, tokens[1], tokens[2]);
				g = th.getMatchGraph(m, 0.5, true);
				//g = th.refineMatchedSpectrum(g, m);
				//scorer.computeScore(g, true, true);
				double score = scorer.compare(th, m);
				System.out.println("score is: " + score);
				currentLine = bf.readLine();
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	
	public static void testTheoreticalSingleLinkedPeptide(){
		Peptide p1 = new Peptide("EGVSKDDAEALKK.2");
		Peptide p2 = new Peptide("MATVSMRDMLKAGVHFGHQTR.2");
		int linkedCharge = 4;
		p1.setPos(new int[]{5});
		p2.setPos(new int[]{11});
		p1.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(p2.getPeptide())+Mass.DSPLINKER_MASS});
		System.out.println("ptm1 is: " + p1.getPtmmasses()[0]);
		p2.setPtmmasses(new double[]{PeptideMassAnalysis.computeMolecularMass(p1.getPeptide())+Mass.DSPLINKER_MASS});	
		System.out.println("ptm2 is: " + p2.getPtmmasses()[0]);
		//TheoreticalSpectrum t1 = new TheoreticalSpectrum(p1, linkedCharge, prefixIons, suffixIons);
		//TheoreticalSpectrum t2 = new TheoreticalSpectrum(p2, linkedCharge, prefixIons, suffixIons);
		TheoreticalSpectrum linkedT12 = new TheoreticalSpectrum(p1, p2, (short)linkedCharge, true);
		List<Peak> linkedT12PList = linkedT12.getPeak();
		for(int i = 0; i < linkedT12PList.size(); i++){
			System.out.println(linkedT12PList.get(i));
		}
		System.out.println("========================================================================");
		p1.setPtmmasses(new double[]{145});
		p2.setPtmmasses(new double[]{145});
		TheoreticalSpectrum linked = new TheoreticalSpectrum(p1, p2, (short)4);
		List<Peak> linkedPeakList = linked.getPeak();
		for(int i = 0; i < linkedPeakList.size(); i++){
			System.out.println(linkedPeakList.get(i));
		}
	}
	
	public static void testTheoreticalLinkedPeptide(){
		//String mixtureFile = "..\\mixture_linked\\tk090204_WKarzai_Tryp.mgf";
		//String tripletFile ="..\\mixture_linked\\triplet.txt";
		String mixtureFile = "..\\mixture_linked\\spectrums_raw.mgf";
		//String tripletFile ="..\\mixture_linked\\triplet_xquest.txt";
		String tripletFile ="..\\mixture_linked\\triplet_selectedsubset.txt";
		SpectrumLib mix = new SpectrumLib(mixtureFile, "MGF");
		mix.toRelativeIntensity();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); //skip header
			String[] tokens;
			Spectrum s1, s2, m;
			TheoreticalSpectrum th;
			while(currentLine != null){
				System.out.println("line is: " + currentLine);
				tokens = currentLine.split("\t");
//				m.filterPeaks((tokens[1].length() + tokens[2].length()) * 5);
				System.out.println("id is : " + tokens[0]+".raw");
				m = mix.getSpectrumById(tokens[0]+".raw");
				//m = mix.getRandomSpectrum();
				System.out.println(m.peptide);
				//m.windowFilterPeaks(10, 25);
				m.computePeakRank();
				currentLine = bf.readLine();
				System.out.println("peptides are: " + tokens[1] + " & " + tokens[2]);
				//testMultipleLysPosition(tokens[1]+".2", tokens[2]+".2", m);
				testSpecificLysPosition(tokens[1]+".2", tokens[2]+".2", m, Integer.parseInt(tokens[3]), Integer.parseInt(tokens[4]));
				//return;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			System.out.println(ioe.getCause());
		}
	}
	

	public static void testMultipleLysPosition(Spectrum s1, Spectrum s2, Spectrum m){
		testMultipleLysPosition(s1.peptide, s2.peptide, m);
	}
	
	public static void testMultipleLysPosition(String peptide1, String peptide2, Spectrum m){
		Peptide p1 = new Peptide(peptide1);
		Peptide p2 = new Peptide(peptide2);
		List<Integer> positions1 = p1.getLysPositions();
		List<Integer> positions2 = p2.getLysPositions();
		//System.out.println("total combination links are: " + positions1.size()*positions2.size());
		TheoreticalSpectrum th;
		for(int i = 0; i < positions1.size(); i++){
			for(int j = 0; j < positions2.size(); j++){
				p1.createDSPLinkerPTM(new int[]{positions1.get(i).intValue()+1});
				p2.createDSPLinkerPTM(new int[]{positions2.get(j).intValue()+1});
				th = new TheoreticalSpectrum(p1, p2, (short)m.charge);
				System.out.println("done creating theo-spect");
				List[] matchedpeaks = th.matchSpectrum(m, 0.5);
				//matchedpeaks = th.refineMatchedSpectrum(matchedpeaks, m);
				th.explainedPeaks2(matchedpeaks, m, 0.5);
				List[] l1 = new List[] {new ArrayList(), new ArrayList()};
				List[] l2 = new List[] {new ArrayList(), new ArrayList()};
				th.splitPeaks(matchedpeaks, l1, l2);
				th.IonsStat(l1, m);
				th.IonsStat(l2, m);
				th.printIonsStatDetail(l1, m);
				//th.printIonsStatDetailWithCharge(l1, m);
				System.out.println();
				th.printIonsStatDetail(l2, m);
				//th.printIonsStatDetailWithCharge(l2, m);
				int[] matchRanks = th.getMaxRanks(2, l1, m);
				int[] matchRanks2 = th.getMaxRanks(2, l2, m);
				System.out.println("max ions matched ranks are " + matchRanks[0] + ", " + matchRanks[1]);
				System.out.println("max ions matched ranks are " + matchRanks2[0] + ", " + matchRanks2[1]);
			}
		}
		
	}
	
	public static void testSpecificLysPosition(String peptide1, String peptide2, Spectrum m, int position1, int position2){
		Peptide p1 = new Peptide(peptide1);
		Peptide p2 = new Peptide(peptide2);
		List<Integer> positions1 = p1.getLysPositions();
		List<Integer> positions2 = p2.getLysPositions();
		//System.out.println("total combination links are: " + positions1.size()*positions2.size());
		TheoreticalSpectrum th;
				p1.createDSPLinkerPTM(new int[]{position1});
				p2.createDSPLinkerPTM(new int[]{position2});
				th = new TheoreticalSpectrum(p1, p2, (short)m.charge);
				System.out.println("done creating theo-spect");
				List[] matchedpeaks = th.matchSpectrum(m, 0.5);
				//matchedpeaks = th.refineMatchedSpectrum(matchedpeaks, m);
				th.explainedPeaks2(matchedpeaks, m, 0.5);
				List[] l1 = new List[] {new ArrayList(), new ArrayList()};
				List[] l2 = new List[] {new ArrayList(), new ArrayList()};
				th.splitPeaks(matchedpeaks, l1, l2);
				th.IonsStat(l1, m);
				th.IonsStat(l2, m);
				th.printIonsStatDetail(l1, m);
				//th.printIonsStatDetailWithCharge(l1, m);
				System.out.println();
				th.printIonsStatDetail(l2, m);
				//th.printIonsStatDetailWithCharge(l2, m);
				int[] matchRanks = th.getMaxRanks(3, l1, m);
				int[] matchRanks2 = th.getMaxRanks(3, l2, m);
				System.out.println("max ions matched ranks are " + matchRanks[0] + ", " + matchRanks[1]);
				System.out.println("max ions matched ranks are " + matchRanks2[0] + ", " + matchRanks2[1]);
		
	}
	
	public static void testMatchingGraph(){
		String spectrumFile = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
		String annotationFile =".\\mixture_linked\\trps\\result.txt";
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
		System.out.println("Done loading spectrum library");
		lib1.annoateSpectrumFromInspectFile(annotationFile);
		System.out.println("Done loading linked-mixture library");
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel.txt";
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		List spectList = lib1.getAllSpectrums();
		Spectrum s;
		TheoreticalSpectrum th;
		long start = (new GregorianCalendar()).getTimeInMillis();
		int count = 0;
		for(int iter = 0; iter < 50; iter++){
			for(int i = 0; i < spectList.size(); i++){
				s = (Spectrum)spectList.get(i);
				if(!s.peptide.contains("Scan")){
					th = new TheoreticalSpectrum(s.peptide);
					//th.matchSpectrum(s, 0.5);
					th.getMatchGraph(s, 0.5);
					count++;
				}
			}
		}
		System.out.println("matching " + count + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testSpectPlot(){
		//String spectrumLibFile = "..\\mixture_linked\\yeast_data\\klc_010908p_yeast-digest.mzXML";
		String spectrumLibFile = "..\\mixture_linked\\Linked_peptide_library\\sumo_lib\\20101008_Sumo_Library_4349_Bo.mzXML";
		String annotationFile = "..\\mixture_linked\\testAnnotation.txt";
		List<String> lines = FileIOUtils.createListFromFile(annotationFile);
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		for(int i = 0; i < lines.size(); i++){
			String[] tokens = lines.get(i).split("\\s+");
			int scan = Integer.parseInt(tokens[2]);
			Spectrum s = reader.getSpectrum(scan);
			s.windowFilterPeaks(10, 25);
			//int charge = Integer.parseInt(tokens[5]);
			int charge = s.charge;
			if(charge < 0 ){
				continue;
			}
			String peptide = tokens[3];//.replaceAll("[\\+\\.0-9]", "");
			System.out.println("peptide is: " + peptide + " charge is : " + s.charge);
			LinkedPeptide lp = new LinkedPeptide(peptide, charge, 1, 6);
			//LinkedPeptide lp = LinkedPeptide.createLinkedPeptide(peptide, charge);
			String outfile = "..\\mixture_linked\\specplotAnnotation\\spectPlotAnnotated_"+scan+"_"+peptide/*lp*/+".txt";
			TheoreticalSpectrum linkedSpect = new TheoreticalSpectrum(lp, lp.getCharge(), false);
			//linkedSpect.shiftSpectrumPPM(100);
			linkedSpect.printSpectPlotAnnotationFile(s, outfile, 0.5);
		}
	}
	
	public static void main(String[] args){
		//testTheoreticalMasses();|
		//testTheoreticalMixture();
		//testTheoreticalMixture2();
		//testTheoreticalLinkedPeptide();
		//testMatchingGraph();
		//testTheoreticalSingleLinkedPeptide();
		testSpectPlot();
	}

}