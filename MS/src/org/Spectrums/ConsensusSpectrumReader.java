package org.Spectrums;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import org.apache.axis.utils.ArrayUtil;

import IO.MZXMLReader;
import Utils.ArrayUtils;

public class ConsensusSpectrumReader implements Iterator<Spectrum>{
	private MZXMLReader reader;
	private List<Spectrum> specList;
	public int numNeighbors=1;
	public int cycle;
	private Spectrum  current;
	private int MaxCacheSize = 300;
	private SortedMap<Integer, Spectrum> cache; //we stored the neighbor scans so do not have to retrieve them constantly
	public double minInt=20;
	public boolean decoyMode = false;
	public double tolerance = 0.05;
	int minScan;
	int maxScan;
	Map<Integer,Spectrum> map;
	
	public ConsensusSpectrumReader(String spectrumFile){
		this(new MZXMLReader(spectrumFile));
//		int count = reader.getSpectrumCount();
//		SortedMap<Integer, Spectrum> specMap = new TreeMap();
//		for(int i = 0; i < count; i++){
//			Spectrum s = reader.getSpectrum(i);
//			if(s != null){
//				s.windowFilterPeaks2(15, 20);
//				specMap.put(s.scanNumber, s);
//			}
//		}
//		this.cache = specMap;
//		System.out.println("Done loading specMap");
	}
	
	public ConsensusSpectrumReader(MZXMLReader reader){
		this.reader = reader;
		this.specList = new ArrayList<Spectrum>();
		this.cache = new TreeMap<Integer, Spectrum>();
		this.minScan = 0;
		this.maxScan = this.reader.getSpectrumCount();
	}
	
	public Spectrum getSpectrum(int scan){
		Spectrum currentSpec = this.reader.getSpectrum(scan);
		getNeighborScans(currentSpec);
		//currentSpec.filterPeaksByIntensity(this.minInt);
		//currentSpec.filterPeaks(1000);
		currentSpec.windowFilterPeaks(15, 25);
		computeConsensus(currentSpec, this.specList);
		return currentSpec;
	}
	
	private Spectrum getConsensusSpectrum(){
		this.current = this.reader.next();
		getNeighborScans(this.current);
		//this.current.filterPeaksByIntensity(this.minInt);
		//this.current.filterPeaks(1000);
		computeConsensus(this.current, this.specList);
		return current;
		
	}
	
	private void getNeighborScans(Spectrum s){
		this.specList.clear();
		this.specList = getNeighborScans(s.scanNumber, this.numNeighbors, this.numNeighbors, this.cycle);
	}
	
	public List<Spectrum> getNeighborScans(int scanNum, int leftSpan, int rightSpan){
		return getNeighborScans(scanNum, leftSpan, rightSpan, this.cycle);
	}
	
	public List<Spectrum> getNeighborScans(int scanNum, int leftSpan, int rightSpan, int cycle){
		//System.out.println("Cycle: " + this.cycle);
		//int Cycle = 35;
		List<Spectrum> specList = new ArrayList<Spectrum>();
		//System.out.println("current Scan: " + currentScan);
		this.specList.clear();
		for(int i = -1*leftSpan; i <= rightSpan; i++){
			int scan = i*cycle + scanNum;
			//System.out.println("Getting scan " + scan);
			if(scan > minScan && scan < maxScan){
				if(this.cache.containsKey(scan)){
					specList.add(this.cache.get(scan));
				}else{
					//System.out.println("cannot find it");
					Spectrum s = this.reader.getSpectrum(scan);
					//s.mergePeaks(s, 0.032);
					if(s!=null){
						s.windowFilterPeaks2(15, 20);		
						s.mergePeaks(s, tolerance);
						//System.out.println("tolerance: " + tolerance);
						//System.out.println(s);
						this.cache.put(s.scanNumber, s);
						if(this.cache.keySet().size() > this.MaxCacheSize){
							this.cache.remove(this.cache.firstKey());
						}
						specList.add(s);
					}
				}
			}
		}
		//System.out.println("neighbor size: " + specList.size());
		return specList;
	}
	
	public double[][] getProjections(Spectrum s, List<Spectrum> specList, double tolerance){
		double[][] projections = new double[s.getPeak().size()][specList.size()];
		for(int i = 0; i < specList.size(); i++){
			Spectrum neigh = specList.get(i);
			//System.out.println(neigh);
			double[] projection = neigh.projectArray(s, tolerance);
			for(int j = 0; j < projection.length; j++){
				projections[j][i] = projection[j];
			}
		}
		return projections;
	}
	
	
	public double getProjectCosine(Spectrum s, Spectrum s1){
		double[] sarry = ArrayUtils.getArray(s);
		double[] tarry = s1.projectArray(s, 0.05);
		ArrayUtils.normalize(sarry);
		ArrayUtils.normalize(tarry);
		return ArrayUtils.dotProd(sarry, tarry);
	}
	
	public double[] getProjectCosine(Spectrum s, List<Spectrum> specList, double tolerance){
		boolean DEBUG = false;
		//finding the middle/'apex' spectrum
		int targetInd = (int)Math.floor((specList.size()/2.0));
		
		//finding precursor info to correlate
		Spectrum precur = new Spectrum();
		ArrayList<Peak> pList = new ArrayList<Peak>();
		Peptide p = new Peptide(s.peptide + "." + s.charge);
		//if(Math.abs(p.getParentmass() - s.parentMass) > 0.05){
		//	System.out.println("precursor not exactly correct\t" + s.spectrumName + "\t" + s.peptide + "\t" + s.parentMass + "\t" + p.getParentmass());
		//}
		//s.parentMass = p.getParentmass();
		pList.add((new Peak(s.parentMass, 10000)));
		precur.setPeaks(pList);
		int MS1Scan = 0;//this.reader.getPrevScan(specList.get(targetInd).scanNumber, 1);
		int K = (specList.size()-1)/2;
		//System.out.println("width: " + K);
		List<Spectrum> surveyScanList = this.getNeighborScans(MS1Scan, 5, 5, this.cycle);

		
		double[][] precursorProjection1 = getProjections(precur, surveyScanList, tolerance);
		double[][] projections = getProjections(s, specList, tolerance);
		double[][] rawprojections = getProjections(s, specList, tolerance);  //raw intensity can get raw intensity without being normalize 
		if(DEBUG){
			System.out.println(s);
		}
		double[][] precursorProjection2 = getProjections(precur, specList, tolerance);
		ArrayUtils.normalize(precursorProjection1[0]);
		ArrayUtils.normalize(precursorProjection2[0]);

		
		//System.out.println("projections-dim: " + projections.length + "\t" + projections[0].length);
		Peak basePeak = s.getPeak().get(0);
		int baseInd = 0;
		SortedMap<Double, Integer> intensityMap = new TreeMap<Double, Integer>();
		//find base/reference peaks for Elution profile correlation
		for(int i = 0; i < s.getPeak().size(); i++){
			Peak next = s.getPeak().get(i);
			//if(projections[i][K] > 0){     //we only compute time-corr for matched peaks
			intensityMap.put(-1*next.getIntensity(), i);
			//}
		}
		baseInd = intensityMap.get(intensityMap.firstKey());
		basePeak = s.getPeak().get(baseInd);
		//transform the spectrum to array
		double[] sarry = ArrayUtils.getArray(s);
		ArrayUtils.normalize(sarry);

		
		//System.out.println("targetInd: " + targetInd);
		double[] tarry = new double[s.getPeak().size()];
		for(int i = 0; i < s.getPeak().size(); i++){
			tarry[i] = projections[i][targetInd];
		}
		ArrayUtils.sqrt(tarry);
		ArrayUtils.normalize(tarry);
		
		//extracting only for top peaks for some statistics
		List<Integer> topIndex = new ArrayList<Integer>();
		int top = 0;
		for(int i = 0; i < sarry.length; i++){
			if(sarry[i] > 0){
				top++;
			}
		}
		//System.out.println(s.peptide + "\tTop:\t" + top);
		if(top < 5){
			top = 5;
		}
		top = 10;
		
		Iterator<Double> keys = intensityMap.keySet().iterator(); 
		for(int i = 0; i < top && keys.hasNext(); i++){
			topIndex.add(intensityMap.get(keys.next()));
		}
		
		
		double[][] tarry2 = new double[projections[0].length][top];
		double[] sarry2 = new double[top];
		for(int i = 0; i < topIndex.size(); i++){
			for(int j = 0; j < projections[0].length; j++){
				tarry2[j][i] = projections[topIndex.get(i)][j];
			}
			sarry2[i] = sarry[topIndex.get(i)];
		}
		//normalize array
		ArrayUtils.normalize(sarry2);  //sarry2 already sqrt above
		for(int j = 0; j < projections[0].length; j++){
			ArrayUtils.sqrt(tarry2[j]);
			ArrayUtils.normalize(tarry2[j]);
		}
		if(DEBUG && ArrayUtils.dotProd(sarry, tarry) 
				!= ArrayUtils.dotProd(sarry2, tarry2[targetInd])){
			System.out.println("similarity not matched");
			System.out.println(ArrayUtils.dotProd(sarry, tarry));
			System.out.println(ArrayUtils.getString(sarry));
			System.out.println(ArrayUtils.getString(tarry));
			System.out.println(ArrayUtils.dotProd(sarry2, tarry2[targetInd]));
			System.out.println(ArrayUtils.getString(sarry2));
			System.out.println(ArrayUtils.getString(tarry2[targetInd]));
		}

		//System.out.println("sarry: " + Arrays.toString(sarry));
		//System.out.println("tarry: " + Arrays.toString(tarry));

		//use only for doing pearson correlation
//		double[][] projectionDev = new double[projections.length][projections[0].length];
//		for(int i = 0; i < projections.length; i++){
//			double centeri = ArrayUtils.average(projections[i]);
//			ArrayUtils.offSet(projections[i], -1*centeri, projectionDev[i]);
//			//System.out.println(getString(projectionDev[i]));
//			ArrayUtils.normalize(projectionDev[i]);
//		}
		
		for(int i = 0; i < projections.length; i++){
			//ArrayUtils.sqrt(projections[i]);
			//this.thresholdArray(projections[i], 0);
			//ArrayUtils.normalize(projections[i]);
			if(DEBUG) System.out.println("parry: " 
					+ s.getPeak().get(i).getMass() + "\t" + s.getPeak().get(i).getIntensity() + "\t"
					+ ArrayUtils.getString(projections[i]));
			ArrayUtils.normalize(projections[i]);
		}
		
		double[] corrs = new double[s.getPeak().size()];
		double[] marry  = new double[s.getPeak().size()];
		double[] marry2  = new double[s.getPeak().size()];
		//computing all pairwise similarity
		double sim = 0.0;
		double baseSim = 0.0; //sim to a base/reference peak
		int count = 0;
		int matchCount = 0;
		int dropCount = 0;
		SortedMap<Double,Integer> corrMap = new TreeMap();
		for(int i = 0; i < s.getPeak().size(); i++){
			double localSim = 0;
			for(int j = 0; j < s.getPeak().size(); j++){
				double corr = ArrayUtils.dotProd(projections[i], projections[j]);
				//corr = corr * sarry[i];
				if(corr > -10 && projections[i][K] > 0 && projections[j][K] > 0){
					sim+=corr;
					localSim+=corr;
					count++;
					//System.out.println("baseInd: " + baseInd);
					if(j == baseInd){
						baseSim +=corr;
						corrs[i] = corr;
					}
				}
			}
			corrs[i] = localSim / projections.length;
			if(corrs[i] > 0 && tarry[i] > sarry[i])
				corrMap.put(corrs[i], i);
			s.computePeakRank();
//			System.out.println(s.spectrumName + "\t" + ArrayUtils.dotProd(sarry, tarry) + "\t" + s.peptide + "\t" 
//				+ s.getPeak().get(i).getIntensity() + "\t" + rawprojections[i][targetInd] + "\t" 
//				+ sarry[i] +"\t" + tarry[i] + "\t"
//				+ s.getPeak().get(i).getRank() + "\t"
//				+ localSim/s.getPeak().size());
		}
		for(int i = 0; i < sarry.length; i++){
			marry[i] = tarry[i];
			if(sarry[i] > 0 && marry[i] > 0){
				matchCount++;
			}
		}
		
		Iterator<Double> it = corrMap.keySet().iterator();
		//System.out.println(ArrayUtils.getString(sarry));
		//System.out.println(ArrayUtils.getString(tarry));
		//System.out.println(ArrayUtils.getString(marry));
		for(int i = 0; i < 1 && it.hasNext(); i++){
			int ind = corrMap.get(it.next());
			if(corrs[ind] < 0.3){
				//System.out.println("setting " + ind + "\t" + corrs[ind] + "\t" +  sarry[ind]);
				marry[ind] = 0;
				sarry[ind] = 0;
			}
		}
		ArrayUtils.normalize(sarry);
		ArrayUtils.normalize(marry);
		//System.out.println(ArrayUtils.getString(sarry));
		//System.out.println(ArrayUtils.getString(tarry));
		//System.out.println(ArrayUtils.getString(marry));

		baseSim /= s.getPeak().size()-1;
			
		double topSim = 0;
		double topBaseSim = 0;
		double topPreSim1 = 0;
		double topPreSim2 = 0;
		int topCount = 0;		
		for(int i = 0; i < topIndex.size(); i++){
			for(int j = i+1; j < topIndex.size(); j++){
				double corr = ArrayUtils.dotProd(projections[topIndex.get(i)], projections[topIndex.get(j)]);
				//corr = corr * sarry[topIndex.get(i)]*sarry[topIndex.get(j)];
				if(corr > -10){
					topSim+=corr;
					topCount += 1;
				}
			}
			if(i == 0){
				topBaseSim = topSim / (topIndex.size()-1);
			}
			topPreSim1 += ArrayUtils.dotProd(projections[topIndex.get(i)], precursorProjection1[0]);
			topPreSim2 += ArrayUtils.dotProd(projections[topIndex.get(i)], precursorProjection2[0]);
		}
		topSim /= topCount;
		topPreSim1 /= top;
		topPreSim2 /= top;

		//compute peak presence
		double topPresence = 0;
		double width = 0;
		for(int i = 0; i < topIndex.size(); i++){
			double[] presencestat =symPresence(projections, topIndex.get(i));
			double corr = presencestat[0];
			//System.out.println("count: " + presencestat[0] + "\t" + corr);
			width = presencestat[1] > width ? presencestat[1] : width;  
			if(DEBUG) System.out.println("top: " + topIndex.get(i));
			if(DEBUG) System.out.println("presence-corr: " + corr);
			//corr /= projections[i].length;
			topPresence+= corr;
			corrs[i] = corr;
		}
		//topPresence /= top;
		
		//topPresence = 0;
		for(int i = 0; i < tarry2.length; i++){
			//double corr = ArrayUtils.dotProd(sarry2, tarry2[i]);
			//if(DEBUG) System.out.println("corr " + corr);
			//topPresence += corr;
		}		
		//topPresence /= tarry2.length;

		double presenceCount  = 0;
		for(int i = 0; i < projections.length; i++){
			double[] presencestat =symPresence(projections, i);
			double corr = presencestat[0];		
			if(DEBUG) System.out.println("presence: " + corr);
			//corr /= projections[i].length;
			//corrs[i] = corr;
			presenceCount +=corr;
			
			//idealize similarity
			//if(corr > 0){
			if(tarry[i] > 0){
				marry2[i] = sarry[i];//corr*tarry[i];
			}else{
				marry2[i] = 0;
			}
		}
		//presenceCount /= projections.length;

		ArrayUtils.normalize(marry2);
		
		if(DEBUG) System.out.println(s.spectrumName + " sarry: " + ArrayUtils.getString(sarry));
		if(DEBUG) System.out.println(s.spectrumName + " tarry: " + ArrayUtils.getString(tarry));
		if(DEBUG) System.out.println(s.spectrumName + " corrs: " + ArrayUtils.getString(corrs));
//		System.out.println(s.spectrumName + " marry: " + ArrayUtils.getString(marry));
		
	
		return new double[]{ArrayUtils.dotProd(sarry, tarry), 
				ratioScore(sarry, tarry), 
				ArrayUtils.dotProd(sarry, marry),
				ArrayUtils.dotProd(sarry, marry2),
				sim/count, topSim, baseSim, topBaseSim, 
				topPresence, presenceCount,
				matchCount, topPreSim1, topPreSim2};
	}
	
	//compute relative intensity of peaks directly
	public static double ratioScore(double[] sarry, double[] tarry){
		double[] rel1 = new double[sarry.length*(sarry.length-1)/2];
		double[] rel2 = new double[tarry.length*(tarry.length-1)/2];
		double[] ratios = new double[sarry.length];
		int k = 0;
		for(int i = 0; i < sarry.length; i++){
			if(sarry[i] < tarry[i]){
				ratios[i] = sarry[i] / (tarry[i] + 0.0000001);
			}else{
				ratios[i] = tarry[i] / (sarry[i] + 0.0000001);
			}
			k++;
		}
//		for(int i = 0; i < sarry.length; i++){
//			for(int j = i+1; j < sarry.length; j++){
//				if(sarry[i] > sarry[j]){
//					rel1[k] = sarry[j] / (sarry[i] + 0.0000001);
//				}else{
//					rel1[k] = sarry[i] / (sarry[j] + 0.0000001);
//				}
//				k++;
//			}
//		}
//		k=0;
//		for(int i = 0; i < tarry.length; i++){
//			for(int j = i+1; j < tarry.length; j++){
//				if(tarry[i] > tarry[j]){
//					rel2[k] = tarry[j] / (tarry[i] + 0.0000001);
//				}else{
//					rel1[k] = tarry[i] / (tarry[j] + 0.0000001);
//				}
//				k++;
//			}
//		}
		double ratioScore = 0;
		for(int i = 0; i < ratios.length; i++){
//			if(rel1[i] > rel2[i]){
//				ratios[i] = rel2[i] / (rel1[i] + 0.000000001);
//			}else{
//				ratios[i] = rel1[i] / (rel2[i] + 0.0000000001);
//			}
			ratioScore += ratios[i]*sarry[i];
		}
		return ratioScore;		
	}
	
	
	public void thresholdArray(double[] arry, double threshold){
		for(int i = 0; i < arry.length; i++){
			if(arry[i] > threshold){
				arry[i] = 1;
			}else{
				arry[i] = 0;
			}
		}
	}
	
	public double[] symPresence(double[][] mat, int ind){
		double count = 0;
		int minWidth = 2;
		int maxWidth = (int)Math.floor(mat[0].length/2);
		int targetInd = maxWidth;
		int width = 2;
		int j = 1;
		for(j = 1; j <= maxWidth; j++){	
			//if(mat[ind][targetInd-j] == 0 
			//		&& mat[ind][targetInd+j] == 0){
				//break;
			//}
			if(targetInd - j > 0 && mat[ind][targetInd-j] > 0){ 
				count++;
			}
			//System.out.println("size: " + mat.length + "\t" +  mat[0].length);
			//System.out.println("inds: " + ind + "\t" + targetInd + "\t" + j);
			if(targetInd + j < mat[0].length && mat[ind][targetInd+j] > 0){ 
				count++;
			}
		}
		width = j > width ? j : width; 
		if(mat[ind][targetInd] > 0){
			count++;
		}
//		System.out.println("width: " + width);
//		if(count > 2*width+1){
//			System.out.println("count too big\t" + count);
//		}
		return new double[]{count, width};

	}
	
	
	public int consecPresence(double[][] mat, int ind){
		int max = 0;
		int series = 0;
		int maxMiss = 10;
		int miss=maxMiss;
		for(int i = 0; i < mat[ind].length; i++){
			if(mat[ind][i] > 0){
				series++;
			}else if(miss > 0 && series > 0){ //must started the series
				//miss--;
			}else{
				//max = series > max ? series : max;
				//series = 0;
				//miss=maxMiss;
			}
		}
		max = series > max ? series : max;
		//return max;
		return series;
	}

	
	private void computeConsensus(Spectrum s, List<Spectrum> specList){
		s.windowFilterPeaks2(20, 25);
		s.mergePeaks(s, 0.05);
		for(int i = 0; i < specList.size(); i++){
			Spectrum current = specList.get(i);
			//current.filterPeaksByIntensity(this.minInt);
			//current.filterPeaks(1000);
			current.windowFilterPeaks(20, 25);
			current.computeConsensus(s, 0.03);
		}
	}
	
	
	
	@Override
	public boolean hasNext() {
		return this.reader.hasNext();
	}

	@Override
	public Spectrum next() {
		return getConsensusSpectrum();
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}
	
	public static  void testConsensus(){
		String spectrumFile = "../mixture_linked/msdata/gringar/swath_development/14344_UPS1_400fm_Ecolilysate_SWATH_5600.mzXML";
		ConsensusSpectrumReader reader = new ConsensusSpectrumReader(spectrumFile);
		reader.numNeighbors = 2;
		reader.minInt = 0;
		reader.decoyMode = true;
		while(reader.hasNext()){
			Spectrum s = reader.next();
			//s = reader.getSpectrum(11624);
			if(s.scanNumber != 11624){
				//continue;
			}
			//System.out.println("consensus: " + s.scanNumber);
			//s.filterPeaksByIntensity(10);
			//s.windowFilterPeaks2(15, 25);
			s.filterPeaksByRankScore(4);
			//s.windowFilterPeaks2(15, 25);
			//System.out.println(s);
			//return;
		}
	}
	
	public static void main(String[] args){
		testConsensus();
	}

}
