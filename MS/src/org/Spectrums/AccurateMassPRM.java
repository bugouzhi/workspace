package org.Spectrums;

import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;

import IO.MZXMLReader;
import Utils.ArrayUtils;

/**
 * PRM spectrum supporting accurate mass tolerance
 * @author Jian
 *
 */
public class AccurateMassPRM extends PRMSpectrum{

	private SpectrumMap spectMap;
	private double[][][] adjAAScores;  //use mass difference between adjecent aa to model mass accuracy
	private double massTolerance = 0.5;
	private double adjMassTolerance = 0.05;
	private int[] validAAMass;
	private double[] prefixOffsets;
	private double[] sufffixOffsets;
	private int[] aaInds;
	private double[] massErrorBins;
	

	/**
	 * For serializable interface 
	 */
	private static final long serialVersionUID = 12321080L;
	
	
	public AccurateMassPRM(Spectrum s, int charge, SpectrumComparator comp){
		this(s, charge, comp, true);
	}
	
	public AccurateMassPRM(Spectrum s, int charge, SpectrumComparator comp, boolean computeTable){
		super(s, charge, comp);
		this.spectMap = new SpectrumMap(this.spectrum);
		this.adjAAScores = new double[this.scoredSpectrum[0].length][20][3];
		//System.out.println("length of ad matraix: " + this.adjAAScores.length);
		this.prefixOffsets = Mass.getIonsOffset(new String[]{"b", "b-H20", "b-NH3"});
		//System.out.println(Arrays.toString(prefixOffsets));
		this.sufffixOffsets = Mass.getIonsOffset(new String[]{"y", "y-H20", "y-NH3"});
		//System.out.println(Arrays.toString(sufffixOffsets));
		this.validAAMass = Mass.getValidAAMasses(0.9995);
		//System.out.println(Arrays.toString(this.validAAMass));
		this.aaInds = Mass.getAAInds();
		this.massErrorBins = new double[103];
		for(int i = 1; i < 101; i++){
			massErrorBins[i]=i*0.01;
		}
		massErrorBins[101]=1001;
		massErrorBins[102] = 2001;
		//System.out.println(Arrays.toString(this.aaInds));
		if(computeTable) computeAdjAATable();
		//for(int i = 0; i < adjAAScores.length; i++){
		//	System.out.println(Arrays.toString(adjAAScores[i]));
		//}
	}
	
	public void computeAdjAATable(){
		double pm = scaleFactor*(spectrum.parentMass*spectrum.charge-Mass.PROTON_MASS*spectrum.charge-Mass.WATER);
		//System.out.println("pm : " + pm);
		for(int m = 0; m < this.adjAAScores.length; m++){
			for(int a = 0; a < this.validAAMass.length; a++){
				for(int c = 1;  c <= spectrum.charge; c++){
					int aa = this.validAAMass[a];
					double expectDiff = ((double)aa)/c;
					for(int i = 0; i < prefixOffsets.length; i++){
					
						double mass1 = (m+(prefixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						double mass2 = (m-aa + (prefixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						if(mass2 < 0){
							continue;
						}
						SortedMap<Double, Peak> matched1 = this.spectMap.getMatchedPeaks(mass1, this.massTolerance);
						SortedMap<Double, Peak> matched2 = this.spectMap.getMatchedPeaks(mass2, this.massTolerance);
//						if(m == 671){
//							System.out.println("matches: " + mass1 + "\t" + mass2 + "\t" +  aa + "\t" + matched1.size() + "\t" + matched2.size());
//						}
						double massDiff = getMassDiff(matched1.values(), matched2.values(), expectDiff);
						if(massDiff < this.adjMassTolerance && m-aa > 0){
							//System.out.println(mass1 + "\t" + mass2 + "\t" + aa);
							adjAAScores[m][a][0] += 12; 
						}else if(massDiff < 1000){
							adjAAScores[m][a][1] += 3;
						}else{
							adjAAScores[m][a][2] += -1;
						}
						
					}
					for(int i = 0; i < this.sufffixOffsets.length; i++){
						double mass1 = (pm-m+(sufffixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						double mass2 = (pm-m+aa +(sufffixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						if(mass2 > pm){
							continue;
						}
						SortedMap<Double, Peak> matched1 = this.spectMap.getMatchedPeaks(mass1, this.massTolerance);
						SortedMap<Double, Peak> matched2 = this.spectMap.getMatchedPeaks(mass2, this.massTolerance);
//						if(m == 671){
//							System.out.println("matches: " + mass1 + "\t" + mass2 + "\t" +  aa + "\t" + matched1.size() + "\t" + matched2.size());
//						}
						double massDiff = getMassDiff(matched2.values(), matched1.values(), expectDiff);
						if(massDiff < this.adjMassTolerance && m-aa > 0){
							//System.out.println(mass1 + "\t" + mass2 + "\t" + aa);
							adjAAScores[m][a][0] += 12;
						}else if(massDiff < 1000){
							adjAAScores[m][a][1] += 3;
						}else{
							adjAAScores[m][a][2] += -1;
						}
					}
					
				}
			}
		}
		System.out.println("Done creating table");
	}
	
	//get adj AA mass accuracy info
	public int[] getAdjAAStat(Peptide p, int[][] stat){
		double pm = scaleFactor*(p.getParentmass()*p.getCharge()-Mass.PROTON_MASS*spectrum.charge-Mass.WATER);
		int[] massInds = this.getMassIndex(p);
		int[] currStat = new int[3];
		//System.out.println("pm : " + pm);
		for(int m = 1; m < massInds.length; m++){
				for(int c = 1;  c <= 1/*spectrum.charge;*/; c++){
					int aa = massInds[m] - massInds[m-1];
					double expectDiff = ((double)aa)/c;
					for(int i = 0; i < prefixOffsets.length; i++){			
						double mass1 = (massInds[m]+(prefixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						double mass2 = (massInds[m]-aa + (prefixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						if(mass2 < 0){
							continue;
						}
						SortedMap<Double, Peak> matched1 = this.spectMap.getMatchedPeaks(mass1, this.massTolerance);
						SortedMap<Double, Peak> matched2 = this.spectMap.getMatchedPeaks(mass2, this.massTolerance);
//						if(m == 671){
//							System.out.println("matches: " + mass1 + "\t" + mass2 + "\t" +  aa + "\t" + matched1.size() + "\t" + matched2.size());
//						}
						double massDiff = getMassDiff(matched1.values(), matched2.values(), expectDiff);
						int ind = ArrayUtils.getIntervalIndex(massDiff, this.massErrorBins);
						stat[i][ind]++;
						
						if(massDiff < this.adjMassTolerance){
							//System.out.println(mass1 + "\t" + mass2 + "\t" + aa);
							currStat[0] += 1; 
						}else if(massDiff < 1000){
							currStat[1] += 1;//0.7;
						}else{
							currStat[2] += 1;//-0.2;
						}
						//System.out.println("Diff: " + massDiff);
					}
					for(int i = 0; i < this.sufffixOffsets.length; i++){
						double mass1 = (pm-massInds[m]+(sufffixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						double mass2 = (pm-massInds[m]+aa +(sufffixOffsets[i]+Mass.PROTON_MASS*(c-1))*scaleFactor)/c;
						if(mass2 > pm){
							continue;
						}
						SortedMap<Double, Peak> matched1 = this.spectMap.getMatchedPeaks(mass1, this.massTolerance);
						SortedMap<Double, Peak> matched2 = this.spectMap.getMatchedPeaks(mass2, this.massTolerance);
//						if(m == 671){
//							System.out.println("matches: " + mass1 + "\t" + mass2 + "\t" +  aa + "\t" + matched1.size() + "\t" + matched2.size());
//						}
						double massDiff = getMassDiff(matched2.values(), matched1.values(), expectDiff);
						int ind = ArrayUtils.getIntervalIndex(massDiff, this.massErrorBins);
						stat[prefixOffsets.length+i][ind]++;
						
						if(massDiff < this.adjMassTolerance){
							//System.out.println(mass1 + "\t" + mass2 + "\t" + aa);
							currStat[0] += 1; 
						}else if(massDiff < 1000){
							currStat[1] += 1;//0.7;
						}else{
							currStat[2] += 1;//-0.2;
						}
						//System.out.println("Diff: " + massDiff);
					}
				}
			}
			return currStat;
	}
	
	public double[][][] getAdjAAScores() {
		return adjAAScores;
	}


	//scores for each pair of adj aa on the peptide
	public double getAdjAAScore(Peptide p){
		double[] scores = this.getAdjAAScores(p);
		return scores[0] + scores[1] + scores[2];
	}
	
	public double[] getAdjAAScores(Peptide p){
		int[] massInd = getMassIndex(p);
		double adjScore=0, adjScore2 = 0, adjScore3 = 0;
		System.out.println("peptide : " + p);
		for(int i = 0; i < p.getPeptide().length(); i++){
			//System.out.println("massindex: " + massInd[i]);										
			int edgeInd = this.aaInds[p.getPeptide().charAt(i) - 'A'];
			//System.out.println("aa Index: " + massInd[i] + "\t" +  edgeInd);
			adjScore += this.adjAAScores[massInd[i]][edgeInd][0];
			adjScore2 += this.adjAAScores[massInd[i]][edgeInd][1];
			adjScore3 += this.adjAAScores[massInd[i]][edgeInd][2];
		}
		return new double[]{adjScore, adjScore2, adjScore3};
	}
	
	private double getMassDiff(Collection<Peak> peaks1, Collection<Peak> peaks2, double offset){
		double massDiff = 1000;
		for(Iterator<Peak> it = peaks1.iterator(); it.hasNext();){
			Peak p = it.next();
			for(Iterator<Peak> it2 = peaks2.iterator(); it2.hasNext();){
				Peak p2 = it2.next(); 
				double diff = Math.abs(p.getMass() - offset - p2.getMass() );
				massDiff = diff < massDiff ? diff : massDiff;
			}
		}
		if(peaks1.size() == 0 && peaks2.size() == 0){
			return 2000;
		}
		return massDiff;
	}
	
	public static void testGetAdjAAStat(){
		String libFile = "../mixture_linked/human_heck_1pepFDR_msgfdb.mgf";
		String result= "../mixture_linked/mixdb_cid_hiaccur_training_decoymatchest.txt";
		String spectrumFile = "../mixture_linked/msdata/Fedor_Mixture_Cid_Hiaccuracy/LOVO_1_1.mzXML";
		//String libFile = "../mixture_linked/msdata/Fedor_Mixture_Cid_Hiaccuracy/Deconvolut/LOVO_1_1_decon.mgf";
		SpectrumLibMap lib = new SpectrumLibMap(libFile, "MGF");
		List<String> results = Utils.FileIOUtils.createListFromFile(result);
		RankBaseScoreLearner pComp = RankBaseScoreLearner.loadComparator("../mixture_linked/Cid_HiAccuracy_model_z_win12_25.o");
		//SpectrumComparator comp = SpectrumUtil.getRankBaseScorer(lib);
		SpectrumComparator comp = new SimpleProbabilisticScorer(pComp);
		
		int[][] adjAAStat = new int[6][103];
		//for(Iterator it = results.iterator(); it.hasNext();){
		for(Iterator it = lib.getAllSpectrums().iterator(); it.hasNext();){
			Spectrum s = (Spectrum)it.next();
			//String line = (String)it.next();
			//String[] tokens = line.split("\\t");
			//if(tokens[1].contains("5547") || line.contains("#") || Double.parseDouble(tokens[11]) < 0.000000015){
			//	continue;
			//}
			//Spectrum s = lib.getSpecByScan(Integer.parseInt(tokens[1]));
			//s.peptide = tokens[2];//.substring(2, tokens[7].length()-2);
			//s.charge = Integer.parseInt(tokens[5]);
			//s.protein = tokens[3];
			s.removePrecursors(0.5);
			s.windowFilterPeaks(15, 25); 
			System.out.println(s.peptide + "\t" + s.charge);
			DecoySpectrumGenerator dg = new DecoySpectrumGenerator();
			if(s.peptide.contains("-") || s.peptide.contains("+")){
				continue;
			}
			//Peptide p = new Peptide(s.peptide, s.charge);
			Peptide p = new Peptide(dg.shuffle(s.peptide), s.charge);
			AccurateMassPRM prmHiAcc = new AccurateMassPRM(s, s.charge, comp, false);
			int[] stat = prmHiAcc.getAdjAAStat(p, adjAAStat);
			System.out.println(s.peptide  + "\t" + Arrays.toString(stat));
			TheoreticalSpectrum th = new TheoreticalSpectrum(p, new String[]{"b", "b-H20", "b-NH3"}, new String[]{"y",  "y-H20", "y-NH3"});
			//th.analyzeAnnotation(s, p.getPeptide(), 0.05);
		}
		for(int i = 0; i < adjAAStat.length; i++){
			System.out.println(Arrays.toString(adjAAStat[i]));	
		}
	}
	
	public static void testCreateAdjTable(){
		String result= "../mixture_linked/mixdb_cid_hiaccur_training_decoymatchest_0.5tol.txt";
		String spectrumFile = "../mixture_linked/msdata/UPS_Ecoli/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		String libFile = "../mixture_linked\\msdata\\Training\\MSGFDB_Tryp_7.mgf";
		List<String> results = Utils.FileIOUtils.createListFromFile(result);
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		SpectrumLibMap lib = new SpectrumLibMap(libFile, "MGF");
		RankBaseScoreLearner pComp = RankBaseScoreLearner.loadComparator("../mixture_linked/yeast_single_model_realannotated_win10_25.o");
		//SpectrumComparator comp = SpectrumUtil.getRankBaseScorer(lib);
		SpectrumComparator comp = new SimpleProbabilisticScorer(pComp);
		((SimpleProbabilisticScorer)comp).matchTolerance =0.5;
		((SimpleProbabilisticScorer)comp).setMinMatchedPeak(0);
		
		SpectrumComparator comp2 = new SimpleProbabilisticScorer(pComp);
		((SimpleProbabilisticScorer)comp2).matchTolerance =0.05;
		((SimpleProbabilisticScorer)comp2).setMinMatchedPeak(0);

		int counter = 0;
		Iterator<String> it = results.iterator();
		//Iterator it = lib.getAllSpectrums().iterator();
		it.next(); //skippng header
		while(it.hasNext()){ 
			String line = (String)it.next();
			String[] tokens = line.split("\\t");
			if(tokens.length < 7){// || Integer.parseInt(tokens[1]) != 41295 || tokens[3].equals("ETD")){
				continue;
			}
			Spectrum s = lib.getSpecByScan(Integer.parseInt(tokens[1]));
			//Spectrum s = (Spectrum)it.next();

			//Spectrum s = reader.getSpectrum(Integer.parseInt(tokens[1]));
			//Spectrum s = (Spectrum)it.next();
			//s.peptide = tokens[2];//tokens[7].substring(2, tokens[7].length()-2);
			//s.charge = Integer.parseInt(tokens[5]);
			s.protein = "RPOTEIN";
			s.removePrecursors(0.5);
			s.windowFilterPeaks(12, 25);
			//s.filterPeaks(150);
			s.computePeakRank();
			Peptide p = new Peptide(s.peptide, 
					s.charge);//tokens[5].substring(tokens[5].length()-1)));
			s.parentMass = p.getParentmass();
			
			System.out.println(s.spectrumName + "\t" + s.scanNumber + "\t" + s.peptide + "\t" + s.parentMass + "\t" + s.charge + "\t" + tokens[2] +"\t" + tokens[6] +"\t" + tokens[3]);
			AccurateMassPRM prmHiAcc = new AccurateMassPRM(s, s.charge, comp);
			double prmScore = prmHiAcc.getScore(p);
			double[] adjScore = prmHiAcc.getAdjAAScores(p);
			s.peptide = p.toString();
			TheoreticalSpectrum th = new TheoreticalSpectrum(p, new String[]{"b", "b(iso)", "b-H20", "b-NH3"}, new String[]{"y", "y(iso)", "y-H20", "y-NH3"});
			SimpleMatchingGraph graph = th.getMatchGraph(s, 0.05);
			SimpleMatchingGraph graph2 = th.getMatchGraph(s, 0.5);
		
			//s.scaleMass(0.9995);
			//th.scaleMass(0.9995);
			//th.analyzeAnnotation(s, p.getPeptide(), 0.5);
			System.out.println(s.spectrumName + "\t" + s.peptide + "\t" + s.charge + "\t" + s.protein + "\tadj-score:\t" 
						+ graph.getVerticeWithEdges(SimpleMatchingGraph.Observed, 1).size() + "\t" 
						+ graph2.getVerticeWithEdges(SimpleMatchingGraph.Observed, 1).size() + "\t"
						+ adjScore[0]  + "\t" + adjScore[1] + "\t" +  adjScore[2] + "\t" + prmScore + "\t"
						+ comp.compare(th, s) +"\t" +  comp2.compare(th, s));  
						//+ "\t" + tokens[10] + "\t" + tokens[11]);
			
			//DecoySpectrumGenerator d = new DecoySpectrumGenerator();
			//String decoy =d.shuffle(p.getPeptide());
			//Peptide dp = new Peptide(decoy, s.charge);
			String decoy = tokens[2];
			Peptide dp = new Peptide(decoy, Integer.parseInt(tokens[5]));
			s.parentMass = dp.getParentmass();
			s.charge = dp.getCharge();
			s.peptide = dp.toString();
			s.protein = tokens[3];
//			
			prmHiAcc = new AccurateMassPRM(s, s.charge, comp);
			double prmScore2 = prmHiAcc.getScore(dp);
			double[] adjScore2 = prmHiAcc.getAdjAAScores(dp);
			TheoreticalSpectrum th2 = new TheoreticalSpectrum(dp, new String[]{"b", "b(iso)", "b-H20", "b-NH3"}, new String[]{"y", "y(iso)", "y-H20", "y-NH3"});
//					
//			
			graph = th2.getMatchGraph(s, 0.05);
			graph2 = th2.getMatchGraph(s, 0.5);

			System.out.println(s.spectrumName + "\t" + decoy + "\t" +  s.charge + "\t" + s.protein + "\t"+ "\tadj-score:\t" 
					+ graph.getVerticeWithEdges(SimpleMatchingGraph.Theoretical, 1).size() + "\t" 
					+ graph2.getVerticeWithEdges(SimpleMatchingGraph.Observed, 1).size() +"\t"
					+ adjScore2[0]  + "\t" + adjScore2[1] + "\t" +  adjScore2[2] + "\t" + prmScore2 + "\t"
					+ comp.compare(th2, s) +"\t" +  comp2.compare(th2, s));  			
		}
		
	}
	
	public static void main(String[] args){
		testCreateAdjTable();
		//testGetAdjAAStat();
	}
	

	
}
