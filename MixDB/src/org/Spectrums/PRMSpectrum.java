package org.Spectrums;

import java.util.Iterator;
import java.util.List;
import java.util.Vector;

/**
 * Compute PRM Spectrum
 * @author Jian Wang
 *
 */
public class PRMSpectrum extends TheoreticalSpectrum{
	public static double MinMass = 0;
	public static double MaxMass = 2000;
	public static double interval = 1.0;
	public static int minCharge = 1;
	public static int maxCharge = 4;
	public static double tolerance = 10.0;
	
	//dummy variable to enable spectrum scoring
	private static int[] ptmPos = new int[0];
	private static double[] ptmMass = new double[0];
	private static Peptide shortPeptide = new Peptide("KKKKKKKKKKK", 1);
	private static Peptide longPeptide = new Peptide("KKKKKKKKKKKKKKKKKKKKKKK", 1);
	
	double[][] scoredSpectrum;
	Spectrum spectrum;
	SpectrumComparator comp;
	int charge = 0;
	double resolution = 1.0;
	
	public PRMSpectrum(Spectrum s, SpectrumComparator comp){
		this(s, s.charge, comp,1.0);
	}
	
	public PRMSpectrum(Spectrum s, int charge, SpectrumComparator comp){
		this(s, charge, comp,1.0);
	}
	
	
	public PRMSpectrum(Spectrum s, int charge, SpectrumComparator comp, double resolution){
		this.spectrum = new Spectrum(s);
		//this.spectrum = s;
		this.spectrum.scaleMass(0.9995);
		this.spectrum.computePeakRank();
		this.comp = comp;
		this.charge = charge;
		this.resolution = resolution;
		computePRMSpectrum();
	}
	
	protected void computePRMSpectrum(){
		double mass = MinMass;
		double parentMass = (this.spectrum.parentMass*this.charge - this.charge*Mass.PROTON_MASS - Mass.WATER)*0.9995;
		//System.out.println("parentmass is: " + parentMass);
		this.scoredSpectrum = new double[3][(int)Math.ceil((parentMass-MinMass+tolerance)/resolution+1)];
		double interval = resolution;
		int counter=0;
		while(mass <= parentMass+this.tolerance){
			double oldParent = this.spectrum.parentMass;
			double currentScore;
			for(int i = 1; i <= 1; i++){
				parentMass = (this.spectrum.parentMass*this.charge - this.charge*Mass.PROTON_MASS - Mass.WATER)*0.9995;
				double complement = parentMass - mass;
				double[][] basemass = new double[][]{{mass},{complement}};
				TheoreticalSpectrum t = getSpectrum(basemass, this.spectrum, this.charge);
				currentScore = this.comp.compare(t, spectrum);
				//System.out.println(mass + "\t" + currentScore);
				//ystem.out.println(mass + this.comp.compare(t, spectrum));
				this.scoredSpectrum[i][counter] = Math.round(currentScore);
			}
			mass+=interval;
			counter++;
		}
		
	}
	
	public double getResolution() {
		return resolution;
	}

	public void setResolution(double resolution) {
		this.resolution = resolution;
	}

	public TheoreticalSpectrum getSpectrum(double[][] basemass, Spectrum s, int charge){
		TheoreticalSpectrum t = new TheoreticalSpectrum();
		List<Peak> theoPeaks = this.generatePeaks(basemass, this.prefixIons, this.suffixIons, ptmPos, ptmMass, 1, this.spectrum.charge);
		//System.out.println("computed peaks: " + theoPeaks);
		t.setPeaks(theoPeaks);
		t.parentMass= s.parentMass;
		t.charge = charge;
		Peptide p = this.guessPeptideLength(s);
		p.setCharge((short)charge);
		t.setPeptide(theoPeaks, guessPeptideLength(s));
		t.p = p;
		return t;
	}
	
	
	private Peptide guessPeptideLength(Spectrum s){
		if(s.parentMass*s.charge < 1310){
			return this.shortPeptide;
		}else{
			return this.longPeptide;
		}
	}
	
	public double[] getScoredSpectrum(double parentmass){
		int diff = (int)Math.round(parentmass - this.spectrum.parentMass*this.spectrum.charge);
		if(diff < 0){
			diff = -1;
		}
		if(diff > 0){
			diff = 1;
		}
		return this.scoredSpectrum[1];
	}
	
	public static void testPRMSpectrum(){
		//SpectrumLib lib = new SpectrumLib(".\\MSPLib\\Lib\\ecoli.msp", "MSP");
		SpectrumLib lib = new SpectrumLib("../mixture_linked/yeast_annotated_spectra.mgf", "MGF");
		lib.removeModSpectra();
		RankBaseScoreLearner pComp = RankBaseScoreLearner.loadComparator("../mixture_linked/yeast_single_peptide_model.o");
		//SpectrumComparator comp = SpectrumUtil.getRankBaseScorer(lib);
		SpectrumComparator comp = new SimpleProbabilisticScorer(pComp);
		((SimpleProbabilisticScorer)comp).setMinMatchedPeak(0);
		Iterator<Spectrum> reader = lib.getAllSpectrums().iterator();
		int counter = 0;
		while(reader.hasNext()){
			Spectrum s = reader.next();
			s.windowFilterPeaks(5, 25);
			s.computePeakRank();
			if(s.scanNumber != 2559){
				continue;
			}
			System.out.println("query: " + s.parentMass + "\t" + s.charge +"\t" + s.peptide);
			PRMSpectrum prmSpect = new PRMSpectrum(s, comp);
			Peptide p = new Peptide(s.peptide);
			double[][] base = prmSpect.computeBaseMass(p.getPeptide(), p.getPos(), p.getPtmmasses());
			prmSpect.computePRMSpectrum();
			for(int i = 0; i < prmSpect.scoredSpectrum.length; i++){
				//System.out.println(s.spectrumName + "\t" + i + "\t" + prmSpect.scoredSpectrum[i]);
			}
			double totalScore = 0;
			double[] scores = prmSpect.getScoredSpectrum(p.getParentmass()*p.getCharge());
			for(int i = 0; i < base[0].length; i++){
				int massIndex = (int)Math.round((0.9995*base[0][i]));
				System.out.println("scored: " + massIndex +"\t" + scores[massIndex]);
				totalScore += scores[massIndex];
			}
			TheoreticalSpectrum t = new TheoreticalSpectrum(s.peptide);
			System.out.println(s.spectrumName + "\t" + s.parentMass + "\t" + s.peptide + "\tTotal score: " + totalScore + "\toriginal score: " + comp.compare(t, s));
			counter++;
			if(counter > 100){
				//return;
			}
		}
	}
	
	public static void main(String[] args){
		testPRMSpectrum();
	}
	
}
