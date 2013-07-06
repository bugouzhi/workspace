package org.Spectrums;

import java.util.Collections;
import java.util.List;
import java.util.Vector;

public class LinkedTheoreticalSpectrum extends TheoreticalSpectrum{
	//for half-linked or dangle linked peptides
	public LinkedTheoreticalSpectrum(Peptide p1, int linkedCharge){
		 this(p1, linkedCharge, new String[]{"b", "b-H20", "b-NH3"}, new String[]{"y", "y-H20", "y-NH3"});
		// this(p1, linkedCharge, new String[]{"b"}, new String[]{"y"});
	}
	
	public LinkedTheoreticalSpectrum(Peptide p1, int linkedCharge, String[] prefix, String[] suffix){
		this.peptide = p1.getPeptide()+"." + linkedCharge;
		p1.setCharge((short)2);
		List theoPeaks = this.generatePeaks(p1.getPeptide(), prefix, suffix, p1.getPos(), p1.getPtmmasses(), p1.getCharge());
		p1.setCharge((short)1);
		//this.shiftPeaks(theoPeaks, p1.getPos(), p1.getPtmmasses(), p1.getPeptide().length());
		
		this.addChargedPeaks(theoPeaks, p1.getPos(), linkedCharge);
		
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
	public LinkedTheoreticalSpectrum(Peptide p1, Peptide p2, short charge){
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
	
	public LinkedTheoreticalSpectrum(Peptide p1, Peptide p2, short charge, boolean dummy){
		super();
		this.peptide = p1 + " & " + p2;
		
		TheoreticalSpectrum th1 = new LinkedTheoreticalSpectrum(p1, charge, this.prefixIons, this.suffixIons);
		TheoreticalSpectrum th2 = new LinkedTheoreticalSpectrum(p2, charge, this.prefixIons, this.suffixIons);
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
		this.parentMass= (p1.getParentmass() + p2.getParentmass() + Mass.DSPLINKER_MASS + charge*Mass.PROTON_MASS)/charge;
		this.charge = charge;
	}
	
	public LinkedTheoreticalSpectrum(TheoreticalSpectrum t1, TheoreticalSpectrum t2, short charge, boolean dummy){
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

}
