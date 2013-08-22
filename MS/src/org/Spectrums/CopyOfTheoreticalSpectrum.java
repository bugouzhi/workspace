package org.Spectrums;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
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
import org.biojava.bio.proteomics.MassCalc;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.ProteinTools;
import org.jgrapht.*;
import org.jgrapht.graph.*;



/**
 * Theoretical Spectrum of a peptide
 * @author jian wang
 *
 */
public class CopyOfTheoreticalSpectrum extends Spectrum{
	public static int MAX_PROFILE_SPAN = 0;
	public  static String[] prefixIons = {"b", "b(iso)", "b-H20"};//, "b-H20-H20", "b-H20-NH3"};//, "a", "a-H20", "a-NH3"};
	public static  String[] suffixIons = {"y", "y(iso)", "y-H20"};//, "y-H20-H20", "y-H20-NH3"};
	//private static String[] prefixIons = { "b" , "b-H20"};//, "a", "a-H20", "a-NH3"};
	//private static  String[] suffixIons = {"y" , "y-H20"};

	private static double MS2Tolerance = 0.5;
	public  Peptide p;
	
	
	protected CopyOfTheoreticalSpectrum(){
		
	}
	
	public CopyOfTheoreticalSpectrum(String peptide){
		this(new Peptide(peptide));
	}
	public CopyOfTheoreticalSpectrum(String peptide, String[] prefixIons, String[] suffixIons){
		this(new Peptide(peptide), prefixIons, suffixIons);
		Vector<Peak> theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
				p.getPos(), p.getPtmmasses(), p.getCharge());
//		Vector theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
//				null, null, p.getCharge());
//		this.shiftPeaks(theoPeaks, p.getPos(), p.getPtmmasses());
		this.setPeaks(theoPeaks);
		this.peptide = p.getPeptide()+"."+p.getCharge();
		this.parentMass= PeptideMassAnalysis.computeMbyZ(peptide, p.getCharge());
		this.charge = p.getCharge();
		this.setPeptide(theoPeaks, p);
		this.p = p;
	}
	public CopyOfTheoreticalSpectrum(Peptide p){
		super();
		Vector<Peak> theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
				p.getPos(), p.getPtmmasses(), p.getCharge());
//		Vector theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
//				null, null, p.getCharge());
//		this.shiftPeaks(theoPeaks, p.getPos(), p.getPtmmasses());
		this.setPeaks(theoPeaks);
		this.peptide = p.getPeptide()+"."+p.getCharge();
		this.parentMass= PeptideMassAnalysis.computeMbyZ(peptide, p.getCharge());
		this.charge = p.getCharge();
		this.setPeptide(theoPeaks, p);
		this.p = p;
//		if(this.p.getPos().length==0){
//			System.out.println("peptide: " + this.peptide);						
//		}else{
//			System.out.println("peptide: " + this.peptide + " +" + this.p.getPtmmasses()[0]);
//		}
//		for(int i = 0; i < this.getPeaks().size(); i++){
//			System.out.println(this.getPeaks().get(i));
//		}
//		System.out.println();
	
	}
	
	public CopyOfTheoreticalSpectrum(Peptide p, String[] prefixIons, String[] suffixIons){
		super();
		Vector<Peak> theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
				p.getPos(), p.getPtmmasses(), p.getCharge());
//		Vector theoPeaks = this.generatePeaks(p.getPeptide(), prefixIons, suffixIons,
//				null, null, p.getCharge());
//		this.shiftPeaks(theoPeaks, p.getPos(), p.getPtmmasses());
		this.setPeaks(theoPeaks);
		this.peptide = p.getPeptide()+"."+p.getCharge();
		this.parentMass= PeptideMassAnalysis.computeMbyZ(peptide, p.getCharge());
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
	public CopyOfTheoreticalSpectrum(String peptide1, String peptide2){
		super();
		Peptide p1 = new Peptide(peptide1);
		Peptide p2 = new Peptide(peptide2);
		//this.peptide = p1.getPeptide() + " & " + p2.getPeptide();
		this.peptide = p1.toString() + " & " + p2.toString();
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
	
	public CopyOfTheoreticalSpectrum(String peptide1, String peptide2, String[] prefixes, String[] suffixes){
		super();
		Peptide p1 = new Peptide(peptide1);
		Peptide p2 = new Peptide(peptide2);
		//this.peptide = p1.getPeptide() + " & " + p2.getPeptide();
		this.peptide = p1.toString() + " & " + p2.toString();
		short charge = (short)(p1.getCharge() + p2.getCharge());
		double[] ptmmass1 = new double[] {0};
		double[] ptmmass2 = new double[] {0};
		//ptmmass1[0] = 0;
		//ptmmass2[0] = 0;

		Vector<Peak> theoPeaks = this.generatePeaks(p1.getPeptide(), prefixes, suffixes,
				p1.getPos(), ptmmass1, p1.getCharge());
		Vector<Peak> theoPeaks2 = this.generatePeaks(p2.getPeptide(), prefixes, suffixes,
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
	public CopyOfTheoreticalSpectrum(Spectrum s1, Spectrum s2){
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
		
		CopyOfTheoreticalSpectrum th1 = new CopyOfTheoreticalSpectrum(s1.peptide);
		CopyOfTheoreticalSpectrum th2 = new CopyOfTheoreticalSpectrum(s2.peptide);

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
	public CopyOfTheoreticalSpectrum(Peptide p1, int linkedCharge){
		 this(p1, linkedCharge, new String[]{"b", "b-H20", "b-NH3"}, new String[]{"y", "y-H20", "y-NH3"});
	}
	
	public CopyOfTheoreticalSpectrum(Peptide p1, int linkedCharge, String[] prefix, String[] suffix){
		this.peptide = p1.getPeptide()+p1.getCharge();
		//System.out.println("ptm pos is : " + p1.getPos()[0]);
		//System.out.println("ptm mass shift is : " + p1.getPtmmasses()[0]);
		List theoPeaks = this.generatePeaks(p1.getPeptide(), prefix, suffix, p1.getPos(), p1.getPtmmasses(), (short)3);
		
		//this.shiftPeaks(theoPeaks, p1.getPos(), p1.getPtmmasses(), p1.getPeptide().length());
		
		this.addChargedPeaks(theoPeaks, p1.getPos(), linkedCharge);
		
		//this.setLinkedPeptide(theoPeaks, p1, linkedCharge);
		p1.setCharge((short)linkedCharge);
		this.setPeptide(theoPeaks, p1);	
		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		this.parentMass = p1.getParentmass(); 
		//this.charge = p1.getCharge();
		this.charge = linkedCharge;
		this.p = p1;
		this.peptide = p1.getPeptide() + "." + p1.getCharge();
//		for(int i = 0; i < this.getPeaks().size(); i++){
//			System.out.println(this.getPeaks().get(i)+"\n");
//		}
	}
	
	//for linked peptide
	public CopyOfTheoreticalSpectrum(Peptide p1, Peptide p2, short charge){
		super();
		this.peptide = p1.getPeptide() + "&" + p2.getPeptide();
		double[] ptmmass1 = new double[] {PeptideMassAnalysis.computeMbyZ(p2.getPeptide(), 1)
				+ Mass.DSPLINKER_MASS -146 }; //we discount an extra hydrogen here since it is include in weight of 2nd peptide
		double[] ptmmass2 = new double[] {PeptideMassAnalysis.computeMbyZ(p1.getPeptide(), 1)
				+ Mass.DSPLINKER_MASS -146 };
		
		CopyOfTheoreticalSpectrum th1 = new CopyOfTheoreticalSpectrum(p1);
		CopyOfTheoreticalSpectrum th2 = new CopyOfTheoreticalSpectrum(p2);

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
	
	public CopyOfTheoreticalSpectrum(Peptide p1, Peptide p2, short charge, boolean dummy){
		super();
		this.peptide = p1.getPeptide() + " & " + p2.getPeptide();
		
		CopyOfTheoreticalSpectrum th1 = new CopyOfTheoreticalSpectrum(p1, charge, this.prefixIons, this.suffixIons);
		CopyOfTheoreticalSpectrum th2 = new CopyOfTheoreticalSpectrum(p2, charge, this.prefixIons, this.suffixIons);

		List matchedPeaks1 = th1.getPeak();//annotatePeaks(th1, s1);
		List matchedPeaks2 = th2.getPeak();//annotatePeaks(th2, s2);
		Vector<Peak> theoPeaks = new Vector();

		setLinkedPeptide(matchedPeaks1, p1, charge);
		setLinkedPeptide(matchedPeaks2, p2, charge);
		
		theoPeaks.addAll(matchedPeaks1);
		theoPeaks.addAll(matchedPeaks2);

		Collections.sort(theoPeaks, new PeakMassComparator());
		//System.out.println("total peaks created: " + theoPeaks.size());
		this.setPeaks(theoPeaks);
		this.parentMass= (p1.getParentmass() + p2.getParentmass() + Mass.DSPLINKER_MASS + charge*Mass.PROTON_MASS)/charge;
		this.charge = charge;
	}
	
	private void setPeptide(List peaks, Peptide p){
		for(int i = 0; i < peaks.size(); i++){
			((LabelledPeak)peaks.get(i)).setPep(p);
		}
	}
	
	private void setLinkedPeptide(List peaks, Peptide pep, int linkedCharge){
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
		//return this.p.toString();
		//return this.p.getPeptide() + "." + this.p.getCharge();
	}
	
	public List<Peak> annotatePeaks(CopyOfTheoreticalSpectrum th, Spectrum s){
		List[] matchedPeaks = th.matchSpectrum(s, CopyOfTheoreticalSpectrum.MS2Tolerance);
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
//		double[][] base = computeBaseMass(peptide, pos, ptmMass);
//		Vector<Peak> theoPeaks = new Vector();
//		//generating prefix ions
//		double[] ions = new double[base[0].length];
//		double[] doubleChargedIons = new double[base[0].length];
//
//		for(int i = 0; i < prefixIons.length; i++){
//			this.addIonMod(base[0], ions, prefixIons[i], 0);
//			createPeaks(theoPeaks, ions, prefixIons[i], (short)1);
//			for(short c = 2; c <= charge; c++){
//				this.addCharge(ions, doubleChargedIons, c, 0);
//				createPeaks(theoPeaks, doubleChargedIons, prefixIons[i], c);
//			}
//		}
//		
//		//generating suffix ions
//		for(int i = 0; i < suffixIons.length; i++){
//			this.addIonMod(base[1], ions, suffixIons[i], 0);
//			createPeaks(theoPeaks, ions, suffixIons[i], (short)1);
//			for(short c = 2; c <= charge; c++){
//				this.addCharge(ions, doubleChargedIons, c, 0);
//				createPeaks(theoPeaks, doubleChargedIons, suffixIons[i], c);
//			}
//		}
//		Collections.sort(theoPeaks, new PeakMassComparator());
//		System.out.println("total ions created is: " + theoPeaks.size());
		return generatePeaks(peptide, prefixIons, suffixIons, pos, ptmMass, 1, charge);
//		return theoPeaks;
		
	}
	
	public Vector<Peak> generatePeaks(String peptide, String[] prefixIons, String[] suffixIons,
			int[] pos, double[] ptmMass, int minCharge, int maxCharge){
		double[][] base = computeBaseMass(peptide, pos, ptmMass);
		Vector<Peak> theoPeaks = new Vector();
		//generating prefix ions
		double[] ions = new double[base[0].length];
		double[] chargedIons = new double[base[0].length];

		for(int i = 0; i < prefixIons.length; i++){
			this.addIonMod(base[0], ions, prefixIons[i], 0);
			for(short c = (short)minCharge; c <= maxCharge; c++){
				this.addCharge(ions, chargedIons, c, 0);
				createPeaks(theoPeaks, chargedIons, prefixIons[i], c);
			}
		}
		
		//generating suffix ions
		for(int i = 0; i < suffixIons.length; i++){
			this.addIonMod(base[1], ions, suffixIons[i], 0);
			for(short c = (short)minCharge; c <= maxCharge; c++){
				this.addCharge(ions, chargedIons, c, 0);
				createPeaks(theoPeaks, chargedIons, suffixIons[i], c);
			}
		}
		Collections.sort(theoPeaks, new PeakMassComparator());
//		System.out.println("total ions created is: " + theoPeaks.size());
		return theoPeaks;
	}
	
	private void createPeaks(List<Peak> pList, double[] masses, String type, short charge){
		LabelledPeak p;
		for(int i = 0; i < masses.length; i++){
			p = new LabelledPeak(masses[i], LabelledPeak.DEFAULT_INTENS, type, (short)(i+1), charge);
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
		suffixMasses[peptide.length()-1] = prefixMasses[peptide.length()-1];
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
			chargedIonMass[startInd+i] = (ionMasses[i] + Mass.PROTON_MASS*(charge-1))/(double)charge;
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
		double total = actual.sumOfPeaks();
		HashMap<Peak, Peak> presented = new HashMap();
		total = total*total;
		double explained = 0.0;
		double intense = 0.0;
		Peak prev = null;
		for(int i = 0; i < matchedPeak[1].size(); i++){
			if(!presented.containsKey(matchedPeak[1].get(i))){ //avoid double counting
				if(matchedPeak[1].get(i).getIntensity() > 0.15){
					intense += matchedPeak[1].get(i).getIntensity();
				}
				explained += matchedPeak[1].get(i).getIntensity();
				presented.put(matchedPeak[1].get(i), matchedPeak[1].get(i));
			}
		}

		//System.out.println("total : " + total + "\texpalined : " + explained);
		System.out.println(actual.spectrumName + "(" + actual.peptide +  " charge: " + actual.charge + ") has explained b/y ions fraction: " + 
				explained/total + " intense b/y: " + intense/total);
		return explained/total;
	}
	
	public double[] IonsStat(List[] peaks, Spectrum s){
		if(peaks[0].size() == 0){
			return new double[]{0,0,0,0};
		}
		double[] b_y_presented = new double[4];
		List<TreeSet> presented = new ArrayList();
		TreeSet<Integer> b = new TreeSet<Integer>();
		TreeSet<Integer> y = new TreeSet<Integer>();

		LabelledPeak l;
		for(int i = 0; i < peaks[0].size(); i++){
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){	
				l = (LabelledPeak)peaks[0].get(i);
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
		System.out.print(s.spectrumName + "(" + ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide() + "): b/y series:" );
		Iterator<Integer> it = b.iterator();
		int prev = 0, current = 0;
		double bSeries = 0;
		double maxBSeries = 0;
		while(it.hasNext()){
			current = it.next().intValue();
			if(current - prev == 1){
				System.out.print("--b"+current);
				bSeries++;
			}else{
				System.out.print("  b"+current);
				bSeries = 0;
			}
			maxBSeries = bSeries > maxBSeries ? bSeries : maxBSeries;
			prev = current;
		}
		System.out.println("\tfraction of total:\t" + b_y_presented[0]);
		it = y.iterator();
		prev = 0;
		current = 0;
		System.out.print(s.spectrumName + "(" + ((LabelledPeak)peaks[0].get(0)).getPep().getPeptide() + ") b/y series: " );
		double ySeries = 0;
		double maxYSeries = 0;
		while(it.hasNext()){
			current = it.next().intValue();
			if(current - prev == 1){
				System.out.print("--y"+current);
				ySeries++;
			}else{
				System.out.print("  y"+current);
				ySeries = 0;
			}
			maxYSeries = ySeries > maxYSeries ? ySeries : maxYSeries;
			prev = current;
		}
		System.out.println( "\tfration of total:\t" + b_y_presented[1]);
		b_y_presented[2] = maxBSeries+1;
		b_y_presented[3] = maxYSeries+1;
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
			if(count <= topN) System.out.print(curr + " (" + p +") ");
			if(count == topN){
				maxRank1 = curr.intValue();
			}
		
		}
		count = 0;
		for(Iterator<Integer> it = matchedPeakRanks[1].keySet().iterator(); it.hasNext(); count++){
			Integer curr = it.next();
			LabelledPeak p = (LabelledPeak)matchedPeakRanks[1].get(curr);
			if(count <= topN) System.out.print(curr  + " (" + p + ") ");
			if(count == topN){
				maxRank2 = curr.intValue();
			}
	
		}
		System.out.println();
		return new int[]{maxRank1, maxRank2};
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
			System.out.println(toBePrint.get(i));
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
							list1.add(pList1.get(i));
							list2.add(pList2.get(j+k));
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
		List<Peak> list1 = new ArrayList<Peak>();
		List<Peak> list2 = new ArrayList<Peak>();
		Spectrum theoretical = this;
		List<Peak> pList1 = theoretical.getPeak(); 
		List<Peak> pList2 = actual.getPeak();
		SimpleMatchingGraph matchingGraph = new SimpleMatchingGraph();
		int i = 0, j = 0, p = 0;
		double m1, m2, diff;
		Peak p1, p2;
		for(int c = 0; c < pList2.size();c++){
			matchingGraph.addVertex(actual.getPeak().get(c), 1);
		}
		
		pList2.add(new Peak(100000, 0)); //we add a end-guard peak to end of the list
		                                 //makes subsquent iteration easier
		for(int c = 0; c < pList1.size(); c++){
			matchingGraph.addVertex(pList1.get(c), 2);
		}
		p1 = pList1.get(0);
		p2 = pList2.get(0);
		for(i = 0; i < pList1.size(); i++){
			p1 = pList1.get(i);
			//System.out.println("matching " + p1 + " starting " + j);
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
					//matchingGraph.addEdge(p2, p1);
					p = j;
					do{
						matchingGraph.addEdge(p2, p1);
						//System.out.println("matching: " + p2 + " ~ " + p1);
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
			if(p.getPep().getPeptide().equals(firstPeptide)){
				l1[0].add(p);
				l1[1].add(p2);
			}else{
				l2[0].add(p);
				l2[1].add(p2);

			}
		}
	}
	
	public static boolean isLinkedPeak(Peptide pep, LabelledPeak p){ 
		boolean prefixLink = p.isPrefixPeak() && p.getPos() >= pep.getPos()[0];
		boolean suffixLink = p.isSuffixPeak() && p.getPos() > pep.getPeptide().length() - pep.getPos()[0];
		return prefixLink || suffixLink;
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
		for(int i = 0; i < pList.size(); i++){
			current = pList.get(i);
			//add charged peaks for all cross-linked species
			for(int j = 0; j < pos.length; j++){
				if(((current.getType().contains("b") || current.getType().contains("a"))&& current.getPos() >= pos[j])
						|| (current.getType().contains("y") && current.getPos() > length-pos[j])){
					if(current.getCharge() == 1){
						count++;
						for(int c = current.getCharge()+1; c <= charge; c++){
							charged = new LabelledPeak((current.getMass() + Mass.PROTON_MASS*(c-1))/(double)c,
									current.getIntensity(), 
									current.getType(),
									current.getPos(),
									(short)c);
							added.add(charged);
						}
					}
					//toBeRemoved.add(current);
				}else{
					break;
				}
			}
		}
		pList.addAll(added);
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
		List[] matching = this.matchSpectrum(s, tolerance);
		double totalFraction = this.explainedPeaks2(matching, s, tolerance);
		//this.splitPeaks(matching, l1, l2);
		double[] b_y_fraction = this.IonsStat(matching, s);
		this.printIonsStatDetailWithCharge(matching, s);
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
		return analyzeMixtureAnnotation(s, p1, p2, 0.3);
	}
	public double[] analyzeMixtureAnnotation(Spectrum s, String p1, String p2, double tolerance){
		List[] matching = this.matchSpectrum(s, tolerance);
		List[] l1 = new List[] {new ArrayList(), new ArrayList()};
		List[] l2 = new List[] {new ArrayList(), new ArrayList()};
		//make sure the order is the same
		
		double totalFraction = this.explainedPeaks2(matching, s, tolerance);
		//System.out.println("peptides " + p1);
		this.splitPeaks(matching, l1, l2, p1.split("\\.")[0], p2.split("\\.")[0]);
		//this.splitPeaks(matching, l1, l2);
		double[] b_y_fraction1 = this.IonsStat(l1, s);
		double[] b_y_fraction2 = this.IonsStat(l2, s);
		double total = s.sumOfPeaks()*s.sumOfPeaks();
		double merror1 = avgError(l1);
		double merror2 = avgError(l2);
		double intensity1 = intensity(l1);
		double intensity2 = intensity(l2);
		//this.printIonsStatDetail(l1, s);
		//System.out.println();
		//this.printIonsStatDetail(l2, s);
		this.printIonsStatDetailWithCharge(l1, s);
		System.out.println();
		this.printIonsStatDetailWithCharge(l2, s);
		double[] stat = new double[13];
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
		stat[11] = intensity1/total;  /// intensity2;
		stat[12] = intensity2/total;
		return stat;
	}
	
	
	public double avgError(List[] peaks){
		LabelledPeak lp;
		Peak p;
		double massDiff = 0.0;
		if(peaks[0].size() == 0){
			return -0.0;
		}
		for(int i = 0; i < peaks[0].size(); i++){
			if(((Peak)peaks[1].get(i)).getIntensity() > 0.00){	
				lp = (LabelledPeak)peaks[0].get(i);
				p = (Peak)peaks[1].get(i);
				massDiff += Math.abs(lp.getMass() - p.getMass());
				//System.out.println("mass error is: " + (lp.getMass() - p.getMass()));
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
			CopyOfTheoreticalSpectrum th;
			Vector v;
			while(currentLine != null){
				th = new CopyOfTheoreticalSpectrum(currentLine);
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
	
	public static void testTheoreticalMasses(){
		String spectrumFile = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
		String annotationFile =".\\mixture_linked\\trps\\result.txt";
		CopyOfTheoreticalSpectrum t;
		String[] prefixIons = {"b"};
		String[] suffixIons = {"y"};
		SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MGF");
		System.out.println("finish loading");
		lib1.annoateSpectrumFromInspectFile(annotationFile);
		lib1.toRelativeIntensity();
		System.out.println("finish annotating");
		Iterator<Spectrum> it = lib1.iterator();
		Spectrum s;
		List[] matchedpeaks;
		while(it.hasNext()){
			s = it.next();
			
			if(!s.peptide.contains("Scan") && /*s.peptide.contains("+") &&*/ (s.charge >= 1)){
				//s.filterPeaks(6*s.peptide.length());
				//s.windowFilterPeaks(6, 25); //inspect parameter
				t = new CopyOfTheoreticalSpectrum(s.peptide);
				matchedpeaks = t.matchSpectrum(s, 0.5);
				t.explainedPeaks2(matchedpeaks, s, 0.5);
				t.IonsStat(matchedpeaks, s);
				t.printIonsStatDetail(matchedpeaks, s);
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
			CopyOfTheoreticalSpectrum th;
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
		//String spectrumFile = ".\\MSPLib\\Lib\\ecoli.msp";
		String mixSpectrumFile = "..\\mixture_compressed\\20070910-A-6mixforMS2_data22.mgf";
		//SpectrumLib lib1 = new SpectrumLib(spectrumFile, "MSP");
		System.out.println("Done loading spectrum library");
		//SpectrumLib mix = lib1.createMix("..\\mixture_linked\\ecoli_mixture2.name", 400, 1, 1, 0, 1, 2000, false);
		SpectrumLib mix = new SpectrumLib(mixSpectrumFile, "MGF");
		System.out.println("Done loading linked-mixture library");
		String tripletFile = "..\\mixture_linked\\mixtureAnnotation.txt";
		String probFile = ".\\data\\IonsScore.txt";
		String noiseModel = ".\\data\\NoiseModel2.txt";
		//lib1.toRelativeIntensity();
		mix.toRelativeIntensity();
		SimplePeakComparator peakscorer = new SimplePeakComparator(probFile, noiseModel);
		SimpleProbabilisticScorer scorer = new SimpleProbabilisticScorer(peakscorer);
		SimpleMatchingGraph g;
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); 
			String[] tokens, tokens2;
			Spectrum s1, s2, m;
			CopyOfTheoreticalSpectrum th;
			while(currentLine != null){
				tokens = currentLine.split("\t");
				//tokens2 = tokens[0].split(" & ");
				m = mix.getSpectra(tokens[0]).get(0);
				m.windowFilterPeaks(10, 25);
				m.computePeakRank();
				th = new CopyOfTheoreticalSpectrum(tokens[1], tokens[2]);
				th.analyzeMixtureAnnotation(m, tokens[1], tokens[2]);
				g = th.getMatchGraph(m, 0.5);
				g = th.refineMatchedSpectrum(g, m);
				scorer.computeScore(g, true, true);
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
		CopyOfTheoreticalSpectrum linkedT12 = new CopyOfTheoreticalSpectrum(p1, p2, (short)linkedCharge, true);
		List<Peak> linkedT12PList = linkedT12.getPeak();
		for(int i = 0; i < linkedT12PList.size(); i++){
			System.out.println(linkedT12PList.get(i));
		}
		System.out.println("========================================================================");
		p1.setPtmmasses(new double[]{145});
		p2.setPtmmasses(new double[]{145});
		CopyOfTheoreticalSpectrum linked = new CopyOfTheoreticalSpectrum(p1, p2, (short)4);
		List<Peak> linkedPeakList = linked.getPeak();
		for(int i = 0; i < linkedPeakList.size(); i++){
			System.out.println(linkedPeakList.get(i));
		}
	}
	
	public static void testTheoreticalLinkedPeptide(){
		//String mixtureFile = "..\\mixture_linked\\tk090204_WKarzai_Tryp.mgf";
		//String tripletFile ="..\\mixture_linked\\triplet.txt";
		String mixtureFile = "..\\mixture_linked\\spectrums_raw.mgf";
		String tripletFile ="..\\mixture_linked\\triplet_xquest.txt";
		//String tripletFile =".\\mixture_linked\\triplet_selectedsubset.txt";
		SpectrumLib mix = new SpectrumLib(mixtureFile, "MGF");
		mix.toRelativeIntensity();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(tripletFile));
			String currentLine = bf.readLine(); //skip header
			String[] tokens;
			Spectrum s1, s2, m;
			CopyOfTheoreticalSpectrum th;
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
		CopyOfTheoreticalSpectrum th;
		for(int i = 0; i < positions1.size(); i++){
			for(int j = 0; j < positions2.size(); j++){
				p1.createDSPLinkerPTM(new int[]{positions1.get(i).intValue()+1});
				p2.createDSPLinkerPTM(new int[]{positions2.get(j).intValue()+1});
				th = new CopyOfTheoreticalSpectrum(p1, p2, (short)m.charge);
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
		CopyOfTheoreticalSpectrum th;
				p1.createDSPLinkerPTM(new int[]{position1});
				p2.createDSPLinkerPTM(new int[]{position2});
				th = new CopyOfTheoreticalSpectrum(p1, p2, (short)m.charge);
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
		CopyOfTheoreticalSpectrum th;
		long start = (new GregorianCalendar()).getTimeInMillis();
		int count = 0;
		for(int iter = 0; iter < 50; iter++){
			for(int i = 0; i < spectList.size(); i++){
				s = (Spectrum)spectList.get(i);
				if(!s.peptide.contains("Scan")){
					th = new CopyOfTheoreticalSpectrum(s.peptide);
					//th.matchSpectrum(s, 0.5);
					th.getMatchGraph(s, 0.5);
					count++;
				}
			}
		}
		System.out.println("matching " + count + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void main(String[] args){
		//testTheoreticalMasses();
		//testTheoreticalMixture();
		//testTheoreticalMixture2();
		testTheoreticalLinkedPeptide();
		//testMatchingGraph();
		//testTheoreticalSingleLinkedPeptide();
	}

}
