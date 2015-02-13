package org.Spectrums;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import IO.MZXMLReader;
import Utils.StringUtils;

/**
 * An extension to MXGF with edge scores that score pair of consecutive
 * aa along a peptides, can handle hi-mass accuracy data
 * @author Jian
 *
 */
public class MXGFe extends MXGF{
	private AccurateMassPRM accuratePRM;
	
	public MXGFe(double maxMass, double maxScore){
		super(maxMass,maxMass, maxScore);
	//	System.out.println(Arrays.toString(this.validAAMass));
	}
	
	
	public AccurateMassPRM getAccuratePRM() {
		return accuratePRM;
	}


	public void setAccuratePRM(AccurateMassPRM accuratePRM) {
		this.accuratePRM = accuratePRM;
	}

	
	public double computeMaxScore(double[] scores, double[][][] adjAAScore){
		double[] maxScore = new double[scores.length];
		//System.out.println("offset: " + scoreOffSet);
		maxScore[0]=0;
		for(int m1 = 1; m1 < scores.length; m1++){
			double max = -1000;
			for(int a = 0; a < this.validAAMass.length; a++){
				int prevMass = m1 - this.validAAMass[a];
				if(prevMass >= 0){
					 double edgeScore = 0;
					 for(int i = 0; m1 < adjAAScore.length && i < adjAAScore[m1][a].length; i++){
						 edgeScore += adjAAScore[m1][a][i];
					 }
					 double current = maxScore[prevMass]+edgeScore;
					 max = current > max ? current : max;
				}
			}
			maxScore[m1]+=max+scores[m1];
		}
		
		double max = -1000;
		for(int i = 0; i < maxScore.length; i ++){
			max = maxScore[i] > max ? maxScore[i] : max;
			//System.out.println("max at : " + i +"\t" + maxScore[i]);
		}
		if(DEBUG) System.out.println("max score at mass: " + maxMass + " is: " + maxScore[maxScore.length-1] + " the max is: " + max);
		return max;
	}
	
	public double computeMinScore(double[] scores, double[][][] adjAAScore){
		double[] minScore = new double[scores.length];
		//System.out.println("offset: " + scoreOffSet);
		minScore[0]=0;
		for(int m1 = 1; m1 < scores.length; m1++){
			double min = 1000;
			for(int a = 0; a < this.validAAMass.length; a++){
				int prevMass = m1 - this.validAAMass[a];
				if(prevMass >= 0){
					double edgeScore = 0;
					 for(int i = 0; m1 < adjAAScore.length && i < adjAAScore[m1][a].length; i++){
						 edgeScore += adjAAScore[m1][a][i];
					 }
					 double current = minScore[prevMass]+edgeScore;
					 min = current < min ? current : min;
				}
			}
			minScore[m1]+=min+scores[m1];
		}
		double min = 1000;
		for(int i = 0; i < minScore.length; i ++){
			min = minScore[i] < min ? minScore[i] : min;
		}
		return min;
	}
	
	
	public void computeEdgeProb(){
		double[][][] adjAAScores = accuratePRM.getAdjAAScores();
		double max = maxMass > maxMass2 ? maxMass : maxMass2;
		this.minScore = this.computeMinScore(new double[this.scoreHi.length], this.accuratePRM.getAdjAAScores())-3;
		this.maxScore = this.computeMaxScore(new double[this.scoreHi.length], this.accuratePRM.getAdjAAScores())+3;
		//System.out.println("min: " + minScore + "\tmax:\t" + maxScore);
		double[][][] DP = new double[(int)Math.ceil(max)][2][(int)Math.ceil(maxScore-minScore+1)];
		int offSet = (int)(-minScore);
		DP[0][0][offSet] = 1;
		int minS = (int)Math.floor(minScore);
		//System.out.println("offset: " + scoreOffSet);
		//System.out.println("maxmass " + maxMass);
		for(int m1 = 0; m1 < maxMass; m1++){
//			if(m1 < adjAAScores.length)
//				System.out.println(Arrays.toString(adjAAScores[m1]));
			for(int s = minS; s < maxScore; s++){
				for(int a = 0; a < this.validAAMass.length; a++){
					int prevMass = m1 - this.validAAMass[a];
					int prevScore = 0;
					if(m1 < adjAAScores.length){
						prevScore = (int)(s - adjAAScores[m1][a][0] - adjAAScores[m1][a][1] - adjAAScores[m1][a][2]);
					}
	 				double prevCount = 0;
					if(prevMass >= 0 && prevScore > minScore && prevScore < maxScore){
						prevScore = prevScore - minS;
						prevCount = DP[prevMass][0][prevScore];
						DP[m1][0][s - minS] += prevCount*0.05;//*this.aaFreq[a];
						if(DEBUG) {
							if(table[m1][0][s - minS] > 0){
								//System.out.println("prob at: " + m1 + "\t" + table[m1][0][s-minS]);
							}
						}
					}
				}
			}
		}
		this.table = DP;
	}
	
	
	public void computeNodeAndEdgeProb(){
		double[][][] adjAAScores = accuratePRM.getAdjAAScores();
		double max = maxMass > maxMass2 ? maxMass : maxMass2;
		this.minScore = this.computeMinScore(this.scoreHi, this.accuratePRM.getAdjAAScores())-3;
		this.maxScore = this.computeMaxScore(this.scoreHi, this.accuratePRM.getAdjAAScores())+3;
		double[][][] DP = new double[(int)Math.ceil(max)][2][(int)Math.ceil(maxScore-minScore+1)];
		int offSet = (int)(-minScore);
		DP[0][0][offSet] = 1;
		int minS = (int)Math.floor(minScore);
		//System.out.println("offset: " + scoreOffSet);
		//System.out.println("maxmass " + maxMass);
		for(int m1 = 0; m1 < maxMass; m1++){
//			if(m1 < adjAAScores.length)
//				System.out.println(Arrays.toString(adjAAScores[m1]));
			for(int s = minS; s < maxScore; s++){
				for(int a = 0; a < this.validAAMass.length; a++){
					int prevMass = m1 - this.validAAMass[a];
					int prevScore = 0;
					if(m1 < adjAAScores.length){
						prevScore = (int)(s - adjAAScores[m1][a][0] - adjAAScores[m1][a][1] - adjAAScores[m1][a][2] - this.scoreHi[m1]);
					}
	 				double prevCount = 0;
					if(prevMass >= 0 && prevScore > minScore && prevScore < maxScore){
						prevScore = prevScore - minS;
						prevCount = DP[prevMass][0][prevScore];
						DP[m1][0][s - minS] += prevCount*0.05;//*this.aaFreq[a];
						if(DEBUG) {
							if(table[m1][0][s - minS] > 0){
								//System.out.println("prob at: " + m1 + "\t" + table[m1][0][s-minS]);
							}
						}
					}
				}
			}
		}
		this.table = DP;
	}
	
	public static double[] computeMSGFe(Spectrum s, String peptide, SpectrumComparator scorer){
		((SimpleProbabilisticScorer)scorer).setMinMatchedPeak(0);
		s.windowFilterPeaks(12, 25);
		s.computePeakRank();
		Peptide p = new Peptide(peptide);
		System.out.println("peptide is: " + p + "\t" + p.getParentmass() +"\t" + p.getCharge() + "\t" + "\tspectrum: " + s.parentMass + "\t" + s.charge);
		s.parentMass = p.getParentmass();
		s.charge = p.getCharge();
		TheoreticalSpectrum t = new TheoreticalSpectrum(p);
		t.analyzeAnnotation(s, peptide, 0.05);
		AccurateMassPRM prmSpect = new AccurateMassPRM(s, p.getCharge(), scorer);
		double score = scorer.compare(t, s);
		double[] scores = prmSpect.getScoredSpectrum(p.getParentmass()*p.getCharge());
		for(int i = 0; i < scores.length; i++){
			scores[i] = scores[i]*5;             //scaling scores to discretize it
		}
		PRMSpectrumComparator sComp = new PRMSpectrumComparator();
		double[][] base = prmSpect.computeBaseMass(p.getPeptide(), p.getPos(), p.getPtmmasses());
		int[] massInds = new int[base[0].length];
		for(int i= 0; i < massInds.length; i++){
			massInds[i] = (int)(Math.round(base[0][i]*0.9995));
		}
		double totalScore = sComp.compare(t, prmSpect);
		double[] adjAAScores = prmSpect.getAdjAAScores(p);
		System.out.println("adj-aa scores: " + Arrays.toString(adjAAScores));
		double adjAAScore = prmSpect.getAdjAAScore(p);
		double pm1 = p.getParentmass() * p.getCharge();
		MXGFe mixgf = new MXGFe(pm1, 200);
		mixgf.setScoreHi(scores);
		mixgf.setScoreLow(scores);
		mixgf.setAccuratePRM(prmSpect);
		mixgf.computeCondPair(massInds);
		double prob = mixgf.getSpecProb(totalScore, (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
		mixgf.computeEdgeProb();
		double adjprob = mixgf.getSpecProb(adjAAScore, (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
		mixgf.computeNodeAndEdgeProb();
		double combinedprob = mixgf.getSpecProb(totalScore+adjAAScore, (p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*p.getCharge()-Mass.WATER));
		return new double[]{score, totalScore, adjAAScore, prob, adjprob, combinedprob};
	}
	
	public static void testComputeMSGFe(String spectrumFile, String annotationFile, String scorerFile){
		List<String> resultLines = Utils.FileIOUtils.createListFromFile(annotationFile);
		Map<Integer, String[]> table = new HashMap<Integer, String[]>();
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		//SpectrumLibMap reader = new SpectrumLibMap(spectrumFile, "MGF");
		PeakComparator comp = RankBaseScoreLearner.loadComparator(scorerFile);
		SpectrumComparator scorer = new SimpleProbabilisticScorer(comp);
		int scanColumn = 1;
		int peptideColumn = 7;
		int chargeInd = 6;
		int counter = 0;
		for(int r = 0;  r < resultLines.size(); r++){
			String result = resultLines.get(r);
			if(result.startsWith("#")){
				continue;
			}
			String[] tokens = result.split("\\t");
			if(Integer.parseInt(tokens[1]) != 24578){
				continue;
			}
//			System.out.print(result);
			Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[scanColumn]));
			//Spectrum s = reader.getSpecByScan(Integer.parseInt(tokens[scanColumn]));
//			String peptide = StringUtils.getPepSeq(tokens[peptideColumn])
//					+ "." + Integer.parseInt(tokens[chargeInd]);
			String peptide = tokens[9];
			double[] stat = computeMSGFe(s, peptide, scorer);
			for(int i = 0; i < tokens.length; i++){
				System.out.print(tokens[i] + "\t");
			}
			for(int i = 0; i < stat.length; i++){
				System.out.print(stat[i] + "\t");
			}
			System.out.println();
			
			TheoreticalSpectrum t = new TheoreticalSpectrum(peptide);
			
			peptide = tokens[10];
			s = s.removeSharePeaks(t, 0.03);
			s.computePeakRank();
			
			stat = computeMSGFe(s, peptide, scorer);
			for(int i = 0; i < tokens.length; i++){
				System.out.print(tokens[i] + "\t");
			}
			for(int i = 0; i < stat.length; i++){
				System.out.print(stat[i] + "\t");
			}
			System.out.println();

		}

	}
	
	public static void main(String[] args){
		//testComputeSingleMXGF();
		//testComputeMSGF();
		args[0] = "../mixture_linked/14344_swathsearch_mixdb_top1_psms.txt";
		args[1] = "../mixture_linked/msdata/UPS_Ecoli/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		//args[1] = "../mixture_linked/msdata/Fedor_Mixture_Cid_Hiaccuracy/Deconvolut/LOVO_1_1_decon.mgf";
		//args[2] = "../mixture_linked/Cid_HiAccuracy_model_z_win12_25.o";
		args[2] = "../mixture_linked/yeast_NIST_lib_singlepep_win12_25.o";
		testComputeMSGFe(args[1], args[0], args[2]);
	}

}
