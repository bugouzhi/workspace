package org.Spectrums;

import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import Utils.FileIOUtils;

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
	public static int LINKEDPREFIXMODE=0;
	public static int LINKEDSUFFIXMODE=1;
	
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
		//computePRMSpectrum();
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
	
	
	protected void computeLinkedPRMSpectrum(int mode){
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
				TheoreticalSpectrum t = getLinkedSpectrum(basemass, this.spectrum, this.charge, mode);
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
	
	
	public TheoreticalSpectrum getLinkedSpectrum(double[][] basemass, Spectrum s, int charge, int mode){
		TheoreticalSpectrum t = new TheoreticalSpectrum();
		List<Peak> theoPeaks = this.generatePeaks(basemass, this.prefixIons, this.suffixIons, ptmPos, ptmMass, 1, this.spectrum.charge);
		//System.out.println("computed peaks: " + theoPeaks);
		t.setPeaks(theoPeaks);
		t.parentMass= s.parentMass;
		t.charge = charge;
		Peptide p = this.guessPeptideLength(s);
		if(mode == PRMSpectrum.LINKEDPREFIXMODE){
			p.setLinkedPos(230); //may need to check how to not-hard code this
		}else{
			p.setLinkedPos(1);  //also check this, but this is more generalizable
		}
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
	
	public static double computeScore(){
		return 0.0;
	}
	
	public static void testLinkedPRMSpectrum(){
		String spectrumFile = "../mixture_linked/linked_peptide_library/sumo_lib/20101008_Sumo_Library_4349_Bo.mzXML";
		String annotationFile = "../mixture_linked/lib_sumo1_sumo_search_with_pyroQ_0.05pm_tolerance_svm.txt";
		String mixtureTraining = "../mixture_linked/lib_sumo_linked_score_model1_win15_25.o";
		SimpleProbabilisticScorer scorer2 = (SimpleProbabilisticScorer)SpectrumUtil.getLinkedSUMOScorer(mixtureTraining);
		
		//SpectrumComparator comp = SpectrumUtil.getRankBaseScorer(lib);
		SpectrumComparator comp = scorer2;
		((SimpleProbabilisticScorer)comp).setMinMatchedPeak(0);
		List<String> results = FileIOUtils.createListFromFile(annotationFile);
		Iterator<String> iter = results.iterator();
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		int counter = 0;
		while(iter.hasNext()){
			String resultLine = iter.next();
			String[] tokens = resultLine.split("\\s+");
			if(!tokens[12].contains("QTGG")){// || !tokens[3].contains("2038")){
				continue;
			}
			Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[3]));
			s.removePrecursors(0.5);
			s.windowFilterPeaks(15, 25);
			s.computePeakRank();
			System.out.println("query: " + s.parentMass + "\t" + s.charge +"\t" + s.peptide);
			String[] peps = new String[2]; 
			Peptide p = new Peptide(tokens[10], 1);
			p.setLinkedPos(p.getPos()[0]);
			Peptide p2 = new Peptide(tokens[12], 1);
			p2.setLinkedPos(6);
			
			String linkedpep = tokens[10]+"--"+tokens[12];
			
			LinkedPeptide lp = new LinkedPeptide(p, p2, s.charge);
			TheoreticalSpectrum sumo = new TheoreticalSpectrum(lp.peptides[1], s.charge);
			//SpectrumUtil.removeAnnotatedPeaks(s, sumo, 0.5);
			
			System.out.println("linkedpep: " + lp);
			
			s.parentMass = lp.getParentmass();
			PRMSpectrum prmSpectPreLink = new PRMSpectrum(s, comp);
			PRMSpectrum prmSpectSuffLink = new PRMSpectrum(s, comp);
			PRMSpectrum prmSpectTag = new PRMSpectrum(s, comp);
			Peptide lp1 = lp.peptides[0];
			Peptide lp2 = lp.peptides[1];
			
			double[][] base = prmSpectPreLink.computeBaseMass(p.getPeptide(), lp1.getPos(), lp1.getPtmmasses());
			double[][] base2 = prmSpectTag.computeBaseMass(p2.getPeptide(), lp2.getPos(), lp2.getPtmmasses());
			
			//setting last massInd to precursor
			base[0][base[0].length-1]=s.parentMass*s.charge-Mass.PROTON_MASS*s.charge-Mass.WATER;
			base2[0][base2[0].length-1]=s.parentMass*s.charge-Mass.PROTON_MASS*s.charge-Mass.WATER;
			
			((LinkedPeptidePeakScoreLearner)((SimpleProbabilisticScorer)comp).comp).peptideMode = 0;
			prmSpectPreLink.computeLinkedPRMSpectrum(PRMSpectrum.LINKEDPREFIXMODE);
			prmSpectSuffLink.computeLinkedPRMSpectrum(PRMSpectrum.LINKEDSUFFIXMODE);
			((LinkedPeptidePeakScoreLearner)((SimpleProbabilisticScorer)comp).comp).peptideMode = 1;
			prmSpectTag.computeLinkedPRMSpectrum(PRMSpectrum.LINKEDPREFIXMODE);
			
			for(int i = 0; i < prmSpectSuffLink.scoredSpectrum.length; i++){
				//System.out.println(s.spectrumName + "\t" + i + "\t" + prmSpect.scoredSpectrum[i]);
			}
			double substrateScore = 0;
			double tagScore = 0;
			double[] scores = prmSpectPreLink.getScoredSpectrum(p.getParentmass()*p.getCharge());
			double[] scores2 = prmSpectSuffLink.getScoredSpectrum(p2.getParentmass()*p2.getCharge());
			double[] scores3 = prmSpectTag.getScoredSpectrum(p.getParentmass()*p.getCharge());
			
			int[] substrateMasses = new int[base[0].length];
			int[] tagMasses = new int[base2[0].length];
			
			for(int i = 0; i < scores2.length; i++){
				//System.out.println(i + "\tscore:\t" + scores2[i]);
			}
			
			System.out.println("linking pos: " + lp1.getLinkedPos());
			//computing score before link
			for(int i = 0; i < lp1.getLinkedPos()-1; i++){
				int massIndex = (int)Math.round((0.9995*base[0][i]));
				substrateMasses[i] = massIndex;
				System.out.println("scored: " + massIndex +"\t" + scores[massIndex]);
				substrateScore += scores[massIndex];
			}
			System.out.println("DONE PREFIX");
			//computing score after link
			for(int i = lp1.getLinkedPos()-1; i < base[0].length; i++){
				int massIndex = (int)Math.round((0.9995*base[0][i]));
				substrateMasses[i] = massIndex;
				System.out.println("scored: " + massIndex +"\t" + scores2[massIndex]);
				substrateScore += scores2[massIndex];
			}
			System.out.println("SCOREING TAG");
			for(int i = 0; i < base2[0].length; i++){
				int massIndex = (int)Math.round((0.9995*base2[0][i]));
				tagMasses[i] = massIndex;
				System.out.println("tag-scored:    " + massIndex +"\t" + scores3[massIndex]);
				tagScore += scores3[massIndex];
			}
			
			((LinkedPeptidePeakScoreLearner)((SimpleProbabilisticScorer)comp).comp).peptideMode = 0;
			TheoreticalSpectrum t1 = new TheoreticalSpectrum(lp.peptides[0], (short)s.charge);
			TheoreticalSpectrum t = new TheoreticalSpectrum(lp.peptides[0], lp.peptides[1], (short)s.charge, true);
			System.out.println(s.spectrumName + "\t" +  s.parentMass + "\t" + lp +  "\t" + lp.getParentmass()
					+ "\tTotal score: " + substrateScore + "\t" + tagScore + "\t" + (substrateScore+tagScore) 
					+ "\toriginal score: " + comp.compare(t1, s) + "\t" + comp.compare(t, s));
						//computing probability
			
			double mass = (s.parentMass*s.charge - s.charge*Mass.PROTON_MASS - Mass.WATER)*0.9995;
			System.out.println("mass is: " + mass);
			MXGF_sumo sumoProb = new MXGF_sumo(mass, scores, scores2, scores3);
			sumoProb.setMasses(substrateMasses, tagMasses);
			sumoProb.initializeSUMO();
			int linkMass = (int)((lp1.getPtmmasses()[0]+128.09496)*0.9995);
			sumoProb.setLinkerMasses(new int[]{linkMass});
			sumoProb.computeProb();
			double[] prob = sumoProb.getSUMOProb(substrateScore, tagScore, mass);
			System.out.println(s.spectrumName + "\t" +  s.parentMass + "\t" + lp +  "\t" + lp.getParentmass()
					+ "\tprobability:\t" + prob[0] + "\t" + prob[1]);
			System.out.println(resultLine +"\t" + prob[0] +"\t" + prob[1]);
			counter++;
			if(counter > 100){
				//return;
			}
		}
	}
	
	public static void main(String[] args){
		//testPRMSpectrum();
		testLinkedPRMSpectrum();
	}
	
}
