package rt;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumLib;
import org.Spectrums.TDAStat;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.collections.bidimap.*;

import sun.swing.UIAction;
import Utils.ArrayUtils;

import java.io.File;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
/**
 * This class align two set of RT between observed data and 
 * library spectra.  It will also automatically determine the RT tolerance
 * window based on the observed data.
 * @author Jian
 *
 */
public class RTAligner implements Serializable{
	//we used two different model to predict rt
	//the global one build a global linear regression model on all rt from identified petides
	//the local one build a local model from closest neighbors (in rt space) of the target peptides
	public static int GLOBAL = 0;
	public static int LOCAL = 1;
	private double[][] data;
	private Set<double[]> dataset;
	private NavigableMap<Double, Double> rtMap;
	private double minCorrelation = 0.95;
	private double numDev = 5.0;
	private double minSpan = 0.8;
	private double maxSpan = 0.99;
	private double minRT = 0;          //minRT--maxRT specify a range of RT we used to perform the linear correlation
	private double maxRT = Double.MAX_VALUE;
	private int predictionMode = GLOBAL; //prediction can be global or local (see above)
	private int numNeighbors=50; //number of neighbors use to build local regression
	
	//residual map
	SortedMap<Double, double[]> residualMap =  new TreeMap<Double, double[]>();
	
	//the master map contain normalized RT for known peptides
	private Map<String, Double> referenceMap;    
	
	//we define a set of reference peptides from the query run with measured RT
	//this set of references are used to align measured RT with those in the master map
	//so peptide RT can be predicted for other peptides in the query run, not in the reference set
    private Map<String, Double> queryRefMap;
    private NavigableMap<Double, String> queryNeighborMap;
	
	//private 
	public boolean DEBUG = false;
	
	SimpleRegression regression;
	
	
	
	/**
	 * Construct aligner directly from RT info, if only global alginment is needed
	 * no need to specify peptide info in the aligner
	 * @param data
	 */
	public RTAligner(double[][] data){
		this(data, Double.MIN_VALUE, Double.MAX_VALUE);
	}
	
	/**
	 * Construct aligner but only using RT info in the range minRT : maxRT
	 * @param data
	 * @param minRT
	 * @param maxRT
	 */
	public RTAligner(double[][] data, double minRT, double maxRT){
		this.minRT = minRT;
		this.maxRT = maxRT;
		initRTAligner(data);
	}
	
	
	private void initRTAligner(double[][] data){
		this.data = data;
		this.regression = new SimpleRegression();
		this.dataset = new HashSet();
		this.rtMap = new TreeMap();
		for(int i = 0; i < data.length; i++){
			this.dataset.add(data[i]);
			this.rtMap.put(data[i][0], data[i][1]);
		}
		//System.out.println("dataset size " + this.dataset.size());
	}
	
	
	/**
	 * Construct a aligner from two spectral library
	 * the library is assumed to contain RT information
	 * @param lib
	 */
	public RTAligner(SpectrumLib lib, SpectrumLib lib2){
		createRTMap(lib, lib2);
	}
	
	
	public RTAligner(Map<String, Double> rtLib, Map<String, Double> rtLib2){
		createRTMap(rtLib, rtLib2);
	}
	
	
	/**
	 * Create  a aligner from a spectral library
	 * if there are multiple spectra for each peptide we just average the RT 
	 * @param lib
	 * @return
	 */
	private void createRTMap(SpectrumLib lib, SpectrumLib lib2){
		HashMap<String, Double> refRTMap = new HashMap<String, Double>(lib.getSpectrumLibrary().keySet().size());
		HashMap<String, Double> queryRTMap = new HashMap<String, Double>(lib2.getSpectrumLibrary().keySet().size());
		//NavigableMap<Double, String> neighborhoodMap = new TreeMap<Double, String>();
		for(Iterator<String> it = lib.getSpectrumLibrary().keySet().iterator(); it.hasNext();){
			String pepKey = it.next();
			List<Spectrum> specList = lib.getSpectra(pepKey);
			double rt = 0.0;
			for(int i = 0; i < specList.size(); i++){
				rt += specList.get(i).rt;
			}
			rt = rt / specList.size();
			refRTMap.put(pepKey, rt);
		}
		
		for(Iterator<String> it = lib2.getSpectrumLibrary().keySet().iterator(); it.hasNext();){
			String pepKey = it.next();
			List<Spectrum> specList = lib2.getSpectra(pepKey);
			double rt = 0.0;
			for(int i = 0; i < specList.size(); i++){
				rt += specList.get(i).rt;
			}
			rt = rt / specList.size();
			queryRTMap.put(pepKey, rt);
		}
		createRTMap(refRTMap, queryRTMap);	
	}
	
	
	/**
	 * Create  a aligner from a pair of RT maps
	 * @param rtLib
	 * @param rtLib2
	 * @return
	 */
	private void createRTMap(Map<String, Double> rtLib, Map<String, Double> rtLib2){
		this.referenceMap = rtLib2;
		//NavigableMap<Double, String> neighborhoodMap = new TreeMap<Double, String>();
		Set<String> intersection = Utils.SetUtils.getIntersect(rtLib.keySet(), rtLib2.keySet());
		this.queryRefMap = new HashMap<String, Double>(intersection.size());
		this.queryNeighborMap = new TreeMap<Double, String>();
		System.out.println("Intersection between reference and query:  " + intersection.size());
		double[][] data = new double[intersection.size()][2];
		int i = 0;
		for(Iterator<String> it = intersection.iterator(); it.hasNext();){
			String pepKey = it.next();
			double rt = rtLib.get(pepKey);
			this.queryRefMap.put(pepKey, rt);
			this.queryNeighborMap.put(rt, pepKey);
			data[i][0]=rt;
			data[i][1] = rtLib2.get(pepKey);
			//System.out.println("paired rt: "  + data[i][0] + " " + data[i][1]);
			i++;
		}
		this.initRTAligner(data);
	}
	
	/**
	 * Compute the regression model for RTs
	 */
	public void computeRegression(){
		for(int i = 0; i < data.length; i++){
			if(data[i][1] >= this.minRT && data[i][1] < this.maxRT){
				//System.out.println(data[i][0] +"\t" + data[i][1]);
				this.regression.addData(data[i][0], data[i][1]);
			}
		}
		//System.out.println("added data " + data.length);
		this.regression.regress();
		int count = 1;
		//System.out.println("Computing regression");
		//DEBUG=true;
		while(count > 0 || this.regression.getR() < this.minCorrelation && this.dataset.size()/this.data.length > this.minSpan){
				count=removeOutliners(numDev);
				if(count == 0){
					removeOutliner(1);
					//count=1;
				}
				if(DEBUG){
					getRegressionStat();
				}
				
		}
		//System.out.println("Done with global regression");
	}
	
	
	public int getPredictionMode() {
		return predictionMode;
	}


	public void setPredictionMode(int predictionMode) {
		this.predictionMode = predictionMode;
	}


	private void getRegressionStat(){
		System.out.println("Correlation is : " + this.regression.getR());
		System.out.println("MSE is : " + this.regression.getMeanSquareError());
		System.out.println("Sum square err is : " + this.regression.getSumSquaredErrors());
		System.out.println("Total Data used : " + this.regression.getN());
		System.out.println("Slope : " + this.regression.getSlope());
		System.out.println("Intercept : " + this.regression.getIntercept());
		getRTDiffInterval();
	}
	
	private SortedMap<Double, double[]> getResidualMap(boolean absolute){
		return getResidualMap(absolute, this.predictionMode);
	}
	
	private SortedMap<Double, double[]> getResidualMap(boolean absolute, int mode){		
		SortedMap<Double, double[]> residualMap =  new TreeMap<Double, double[]>();
		int i = 0;
		for(Iterator<double[]> it = this.dataset.iterator(); it.hasNext();){
			double[] d = it.next();
			double diff =  this.regression.predict(d[0]) - d[1];
			if(DEBUG){
				System.out.println("residueal:\t" + d[0] + "\t" + diff + "\t" + Math.abs(diff));
			}
			if(mode == RTAligner.LOCAL){
				diff =  this.getLocalAlignedRT(d[0], this.numNeighbors) - d[1];
				//System.out.println("finished " + i);
				System.out.println(this.queryNeighborMap.get(d[0]) + "\tresidueal:\t" + d[0] + "\t" + (this.regression.predict(d[0]) - d[1]) +"\t" + diff +"\t"
								+ Math.abs(this.regression.predict(d[0]) - d[1]) +"\t" + Math.abs(diff));
			}
			if(absolute) diff = Math.abs(diff);
			residualMap.put(diff, d);
			i++;
			
		}
		return residualMap;
	}
	
	//note currently remove outliner only apply for the global regression model
	//may consider also doing this for the local regression models to improve performance
	/**
	 * Remove outliner from the regression data
	 * @param num of outliner to remove
	 */
	public void removeOutliner(int num){
		SortedMap<Double, double[]> residualMap = this.getResidualMap(true, RTAligner.GLOBAL);
		for(int i = 0; i < num; i++){
			double[] remove = residualMap.get(residualMap.lastKey());
			if(DEBUG) System.out.println("removing " + remove[0] +"\t" + remove[1]);
			residualMap.remove(residualMap.lastKey());
			this.regression.removeData(remove[0], remove[1]);
			this.dataset.remove(remove);
		}
	}
	
	/**
	 * Remove outliners that is numDev deviation away from median
	 * @param numDev
	 */
	public int removeOutliners(double numDev){
		SortedMap<Double, double[]> residualMap = this.getResidualMap(false, RTAligner.GLOBAL);
		double[] stat = getResidualStat(residualMap);
		if(DEBUG) System.out.println("size of map " + residualMap.size());
		double median = stat[0];
		double inQuantile = stat[2]-stat[1];
		if(DEBUG) System.out.println("median " + median + "\tInnerQ-range: " + inQuantile);
		int count=0;
		double last = residualMap.lastKey();
		double first = residualMap.firstKey();
		if(Math.abs(last-median)/inQuantile > numDev){
				double[] outliner = residualMap.get(last);
				if(DEBUG) System.out.println("removing: " + outliner[0] + "\t" + outliner[1]);
				residualMap.remove(last);
				this.dataset.remove(outliner);
				this.regression.removeData(outliner[0], outliner[1]);
				count++;
		}
		if(Math.abs(median-first)/inQuantile > numDev){
			double[] outliner = residualMap.get(first);
			if(DEBUG) System.out.println("removing: " + outliner[0] + "\t" + outliner[1]);
			residualMap.remove(first);
			this.dataset.remove(outliner);
			this.regression.removeData(outliner[0], outliner[1]);
			count++;
		}
		return count;
	}
	
	//return stats about the inner-quantile range as well as the range that cover maxSpan fraction of data
	private double[] getResidualStat(SortedMap<Double, double[]> residuals){
		int i = 0;
		Double[] diffs = new Double[residuals.size()];
		for(Iterator<Double> it = residuals.keySet().iterator(); it.hasNext();){
			diffs[i++] = it.next();
		}
		
		double median = diffs[(int)(diffs.length/2.0)];
		double inQuantile = diffs[(int)(diffs.length*0.75)] - diffs[(int)(diffs.length*0.25)];
		double minFract = (1 - this.maxSpan)/2.0;
		return new double[]{median, 
				diffs[(int)(diffs.length*0.25)], // inner quantile 
				diffs[(int)(diffs.length*0.75)],
				diffs[(int)(diffs.length*minFract)], //range that cover maxSpan fraction of data
				diffs[(int)(diffs.length*(1-minFract))]};
	}
	
	
	
	public double[] getRTDiffInterval(){
		return getRTDiffInterval(this.maxSpan);
	}
	
	/**
	 * Obtain RT specific interval rather than a global one for all RTs
	 * This means that RT prediction maybe more or less accurate depending
	 * on the number of data point at a particular RT ranges, the prediction
	 * along the beginning or end of the chromatrography also tends to be more
	 * varaible, so RT precition maybe less accurate.  This provide the use of different
	 * expected error window for RT prediciton depending on the locatio along the 
	 * chromatography
	 * @param currRT
	 * @param fract
	 * @return
	 */
	public double[] getRTInterval(double currRT, double fract){
		double[][] localRTs = getKNeighborRTs(currRT, numNeighbors);
		for(int i = 0; i < localRTs.length; i++){
			if(DEBUG) System.out.println(localRTs[i][0] + "\t" + localRTs[i][1]);
		}
		RTAligner localAlign = new RTAligner(localRTs);
		//localAlign.numDev = 3.0;	
		//localAlign.DEBUG=true;
		localAlign.computeRegression();
		//localAlign.regression = this.regression;
		double predicted = localAlign.getAlignedRT(currRT);
		double[] rtDiff=localAlign.getRTDiffInterval(fract);
		//return new double[]{predicted+rtDiff[0], predicted+rtDiff[1], predicted};
		return rtDiff;
	}
	
	/**
	 * This provide a global expected RT errors between observed RT and predicted RT
	 * for the whole run.  The RT error ranges must cover fract of all data points
	 * @param fract
	 * @return
	 */
	public double[] getRTDiffInterval(double fract){
		SortedMap<Double, double[]> residualMap = this.getResidualMap(false);
		double[] stat = getResidualStat(residualMap);
		double minFract = (1 - this.maxSpan)/2.0;
		//System.out.println("#Median\t0.25Q\t0.75Q\t"+this.maxSpan+"Q\t"+(1-this.maxSpan)+"Q\n");
		//System.out.println(Arrays.toString(stat));
		return new double[]{stat[3], stat[4]};
	}
	
	public double getAlignedRT(double rt){
		return getAlignedRT(rt, this.predictionMode);
	}
	
	public double getAlignedRT(double rt, int mode){
		if(mode == RTAligner.GLOBAL){
			return this.regression.predict(rt);
		}else if(mode == RTAligner.LOCAL){
			return getLocalAlignedRT(rt, this.numNeighbors);
		}
		return 0.0;
	}
	
	private double getLocalAlignedRT(double rt, int numNeighbors){
		double[][] localRTs = getKNeighborRTs(rt, numNeighbors);
		//double[][] localRTs = getNeighborRTs(rt-30, rt+30);
		//System.out.println(Arrays.toString(localRTs));
//		SimpleRegression localRegression = new SimpleRegression();
//		localRegression.addData(localRTs);
//		localRegression.regress();
//		if(rt > 82.4 || rt < 82.3){
//			return 0;
//		}
//		for(int i = 0; i < localRTs.length; i++){
//			System.out.println(localRTs[i][0] + "\t" + localRTs[i][1]);
//		}
		
		RTAligner localAlign = new RTAligner(localRTs);
		localAlign.numDev = 3.0;	
		localAlign.computeRegression();
		SimpleRegression localRegression = localAlign.regression;
		//System.out.println("predicted: " + localRegression.predict(rt));
		return localRegression.predict(rt);
//		double avg = 0;
//		double totalweight=0;
//		for(int i = 0; i < localRTs.length; i++){
//			totalweight+=1/(Math.abs(localRTs[i][0]-rt));
//		}
//		//System.out.println("total weight " + totalweight);
//		for(int i = 0; i < localRTs.length; i++){
//			avg+=(localRTs[i][1]);///(Math.abs(localRTs[i][0]-rt))/totalweight);
//		}
//		avg=avg/localRTs.length;
//		return avg;
//		double[] values = new double[localRTs.length];
//		for(int i = 0; i < localRTs.length; i++){
//			values[i]=localRTs[i][1];
//		}
//		double median =  Utils.ArrayUtils.median(values);
//		//System.out.println("median is : " + median);
//		return median;
	}
	
	
	private double[][] getNeighborRTs(double leftRT, double rightRT){
		SortedMap<Double, String> neighbors=this.queryNeighborMap.subMap(leftRT, rightRT);
		double[][] localRTs = new double[neighbors.size()][2];
		int i = 0;
		for(Iterator<Double> it = neighbors.keySet().iterator(); it.hasNext();){
			double rt = it.next();
			String pepKey = this.queryNeighborMap.get(rt);
			localRTs[i][0] = this.queryRefMap.get(pepKey);
			localRTs[i][1] = this.referenceMap.get(pepKey);
			i++;
		}
		return localRTs;
	}
	
	private double[][] getKNeighborRTs(double rt, int K){
		List<String> neighs = new ArrayList<String>();
		int numNeigh = 0;
		double currDiff = 0;
		int direction = 0;
		double[][] localRTs = new double[K][2];
		
		//get one lower and one higher neighbor to start
		Double leftNeigh = this.queryNeighborMap.lowerKey(rt);
		Double rightNeigh = this.queryNeighborMap.higherKey(rt);
		
		if(leftNeigh == null && rightNeigh == null){
			return localRTs;
		}
		
		
		if(leftNeigh == null){
			direction = 1;
		}else if(rightNeigh == null){
			direction = -1;
		}else{		
			direction = rt - leftNeigh < rightNeigh - rt ? -1 : 1; 
		}
		
		if(direction < 0){
			currDiff = rt - leftNeigh;
			neighs.add(this.queryNeighborMap.get(leftNeigh));
		}
		if(direction > 0){
			currDiff = rightNeigh - rt;
			neighs.add(this.queryNeighborMap.get(rightNeigh));
		}
		
		System.out.println("current is : " + rt);	
		while(numNeigh < K){
			if(leftNeigh==null && rightNeigh==null){
	
				return localRTs;
			}
			if(direction < 0){
				//System.out.println("left " + leftNeigh);
				if(leftNeigh == null){
					direction = 1;
					continue;
				}
				double diff = rt - leftNeigh;
				direction = diff > currDiff ? 1 : -1;
				currDiff = rt - leftNeigh;
				String pepKey = this.queryNeighborMap.get(leftNeigh);
				localRTs[numNeigh][0]=leftNeigh;
				localRTs[numNeigh][1]=this.referenceMap.get(pepKey);
				neighs.add(this.queryNeighborMap.get(leftNeigh));
				leftNeigh = this.queryNeighborMap.lowerKey(leftNeigh);
			}else if(direction > 0){
				//System.out.println("right " + rightNeigh);
				if(rightNeigh == null){
					direction = -1;
					continue;
				}
				double diff = rightNeigh - rt;
				direction =  diff > currDiff ? -1 : 1;
				currDiff = rightNeigh - rt;
				String pepKey = this.queryNeighborMap.get(rightNeigh);
				localRTs[numNeigh][0]=rightNeigh;
				localRTs[numNeigh][1]=this.referenceMap.get(pepKey);
				neighs.add(this.queryNeighborMap.get(rightNeigh));
				rightNeigh = this.queryNeighborMap.higherKey(rightNeigh);
			}
			numNeigh++;
			
		}
		//System.out.println("query rt: " + rt);
//		for(int i = 0; i < localRTs.length; i++){
//			System.out.println(Arrays.toString(localRTs[i]));
//		}
//		System.out.println("");
		return localRTs;
	}
	
	
	
	
	public double getMinRT() {
		return minRT;
	}

	public void setMinRT(double minRT) {
		this.minRT = minRT;
	}

	public double getMaxRT() {
		return maxRT;
	}

	public void setMaxRT(double maxRT) {
		this.maxRT = maxRT;
	}

	/**
	 * Run the RT alignment to aligned RT between library and data
	 * @param refResultFile -- contain IDs that use as reference to perform RT fit (i.e. IDs from previous round of search)
	 * @param resultFile -- raw search result file from MSPLIT-DIA
	 * @param resultWithRT -- resultFile with RT updated using the regression calculator
	 * @param RTInd1
	 * @param RTInd2
	 */
	public static void filterWithRTAlign(String refResultFile, String resultFile, String resultWithRT, int RTInd1, int RTInd2, int pepInd, int chargeInd, double minRT, double maxRT){
		RTAligner align = createLocalAlignerFromResult(refResultFile, RTInd1, RTInd2, pepInd, chargeInd, minRT, maxRT);
		align.computeRegression();
		align.setPredictionMode(RTAligner.LOCAL);
		double[] interval = align.getRTDiffInterval();
		System.out.println("RT tolerance interval: " + interval[0] + "\t" + interval[1]);
		List<String> results = Utils.FileIOUtils.createListFromFile(resultFile);
		PrintStream rtOut = Utils.FileIOUtils.getOutStream(resultWithRT);
		for(int i = 0; i < results.size(); i++){
			String[] tokens = results.get(i).split("\\t");
			if(tokens.length < 33 || tokens[0].startsWith("#")){
				if(tokens[0].startsWith("#")){
					rtOut.println(results.get(i)); //print out header
				}
				continue;
			}
			//System.out.println("line is: " + results.get(i));
			try{
				double RT1 = Double.parseDouble(tokens[RTInd1]);
				double RT2 = Double.parseDouble(tokens[RTInd2]);
				double predicted = align.getAlignedRT(RT1);
				//System.out.println("RT: " + RT1 + "\t" + RT2);
				//if(RT1 < 122.2 || RT1 > 122.4){
				//	continue;
				//}
				double[] localinterval = align.getRTInterval(RT1, 0.99);
				System.out.println(RT1 + "\t" + "interval: " + interval[0] + "\t" + interval[1] +"\t" +tokens[8]);
				int pass = 0;
				int localpass = 0;
				if(predicted-RT2 > interval[0] && predicted - RT2 < interval[1] || RT2 < align.minRT || RT2 > align.maxRT){
					rtOut.println(results.get(i));
					pass = 1;
				}
				if(predicted-RT2 > localinterval[0] && predicted - RT2 < localinterval[1] || RT2 < align.minRT || RT2 > align.maxRT){
					//rtOut.println(results.get(i));
					localpass = 1;
				}
				System.out.println(RT1 + "\t" + RT2 + "\t" + predicted + "\tpassing interval:\t" + localinterval[0] + "\t" + localinterval[1] +"\t" +tokens[8] +"\t" + tokens[pepInd] +"\t" + pass + "\t" + localpass +"\t" + tokens[31] + "\t" + localinterval[1]);
			}catch(NumberFormatException e){
				System.err.println("Warining: RT info for the following line is not valid or available");
				System.err.println(results.get(i));
				continue;
			}
		}
		rtOut.flush();
		rtOut.close();
	}
	
	public static RTAligner createAlignerFromResult(String refResultFile, int RTInd1, int RTInd2, int pepInd, int chargeInd, double minRT, double maxRT){
		List<String> filterResults = Utils.FileIOUtils.createListFromFile(refResultFile);
		Map<String, double[]> RTMap = new HashMap<String, double[]>();
		for(int i = 0; i < filterResults.size(); i++){
			//System.out.println("line is: " + results.get(i));
			String[] tokens = filterResults.get(i).split("\\t");
			if(tokens.length < 30 || tokens[0].startsWith("#") || tokens[pepInd].contains("+")){  //mod peptide is tricky we exclude them from RT alignment for now
				continue;
			}
			String key = tokens[pepInd]+"@"+tokens[chargeInd];
			double[] RTs;
			if(RTMap.containsKey(key)){
				RTs = RTMap.get(key);
				RTs[0] += Double.parseDouble(tokens[RTInd1]);
				RTs[1] += Double.parseDouble(tokens[RTInd2]);
				RTs[0] = RTs[0] / 2;
				RTs[1] = RTs[1] / 2;
			}else{
				RTs = new double[2];
				RTs[0] = Double.parseDouble(tokens[RTInd1]);
				RTs[1] = Double.parseDouble(tokens[RTInd2]);
			}
			RTMap.put(key, RTs);
		}
		double[][] data = new double[RTMap.keySet().size()][2];
		int i = 0;
		for(Iterator<double[]> it = RTMap.values().iterator(); it.hasNext();){
			double[] rts = it.next();
			data[i++] = rts;
		}
		RTAligner align = new RTAligner(data, minRT, maxRT);	
		return align;
	}
	
	public static RTAligner createLocalAlignerFromResult(String refResultFile, int RTInd1, int RTInd2, int pepInd, int chargeInd, double minRT, double maxRT){
		List<String> filterResults = Utils.FileIOUtils.createListFromFile(refResultFile);
		Map<String, Double> RTMap = new HashMap<String, Double>();
		Map<String, Double> RTMap2 = new HashMap<String, Double>();
		for(int i = 0; i < filterResults.size(); i++){
			//System.out.println("line is: " + filterResults.get(i));
			String[] tokens = filterResults.get(i).split("\\t");
			if(tokens.length < 30 || tokens[0].startsWith("#") || tokens[pepInd].contains("+")){  //mod peptide is tricky we exclude them from RT alignment for now
				continue;
			}
			String key = tokens[pepInd]+"@"+tokens[chargeInd];
			double[] RTs;
			if(RTMap.containsKey(key)){
				RTs = new double[]{RTMap.get(key), RTMap2.get(key)};
				RTs[0] += Double.parseDouble(tokens[RTInd1]);
				RTs[1] += Double.parseDouble(tokens[RTInd2]);
				RTs[0] = RTs[0] / 2;
				RTs[1] = RTs[1] / 2;
			}else{
				RTs = new double[2];
				RTs[0] = Double.parseDouble(tokens[RTInd1]);
				RTs[1] = Double.parseDouble(tokens[RTInd2]);
			}
			RTMap.put(key, RTs[0]);
			RTMap2.put(key, RTs[1]);
			System.out.println("Reference RT points\t" + RTs[0] + "\t" + RTs[1]);
		}

		RTAligner align = new RTAligner(RTMap, RTMap2);	
		return align;
	}
	
	
	public static void filterMSPLITDIAWithRT(String resultFile, String outFile, double fdr, double minRT, double maxRT){
		String name = Utils.FileIOUtils.stripExtension(outFile);
		String refResultFile = name +"_refResult_raw.txt";
		String rtResultFile = name + "_rtFiltered.txt";
		String filteredrtResultFile = outFile;// name + "_rtfiltered_fdrfiltered.txt";
		TDAStat stat = new TDAStat(resultFile, 1, 4, 6, 8, 32, -1);
		stat.printResultWithFDRInfo(refResultFile, fdr, TDAStat.PRINT_UNIQUE_PEP);
		System.out.println("Done generating reference");
		filterWithRTAlign(refResultFile, resultFile, rtResultFile, 17, 16, 4, 6, minRT, maxRT);
		System.out.println("Done with RT algined");
		stat = new TDAStat(rtResultFile, 1, 4, 6, 8, 32, -1);
		//delete temporary file
		File f1 = new File(refResultFile);
		File f2 = new File(rtResultFile);
		//System.out.println("file exists: " + f1.exists());
		//System.out.println("file exists: " + f2.exists());
		stat.printResultWithFDRInfo(filteredrtResultFile, fdr);
		//f1.delete();
		//f2.delete();
	}
	
	public static void testRTAligner(){
		String resultFile = "..//mixture_linked//SWATH/testRTConstraints//APSWATH_EIF4a2_testRTFit_refResult.txt";
		RTAligner align = createAlignerFromResult(resultFile, 16, 17, 4, 6, 0, Double.MAX_VALUE);
		align.computeRegression();
		double[] interval = align.getRTDiffInterval();
		System.out.println("RT tolerance interval: " + interval[0] + " : " + interval[1]);
	}
	
	public static void testRunRTAlign(){
		String refResultFile = "..//mixture_linked//SWATH/testRTConstraints//APSWATH_EIF4a2_testRTFit_refResult_raw.txt";
		String resultFile = "..//mixture_linked//SWATH/testRTConstraints//APSWATH_EIF4a2_testRTFit_result.txt";
		String rtResultFile = "..//mixture_linked//SWATH/testRTConstraints//APSWATH_EIF4a2_testRTFit_rtFiltered.txt";
		filterWithRTAlign(refResultFile, resultFile, rtResultFile, 17, 16, 4, 6, 0, Double.MAX_VALUE);
	}
	
	public static void testFilterMSPLITDIAwithRT(){
		String resultFile = "..//mixture_linked//SWATH/testRTConstraints//APSWATH_EIF4a2_r1_phlib_swathmsplitout_withlibiRT.txt";
		String outFile = "..//mixture_linked//SWATH/testRTConstraints//AAPSWATH_EIF4a2_r1_phlib_swathmsplitout_withlibiRT_filter.txt";
		filterMSPLITDIAWithRT(resultFile, outFile, 0.01,  Double.MIN_VALUE, Double.MAX_VALUE);
	}
	
	public static void testLocalRTAlign(){
		String resultFile = "..//mixture_linked//SWATH/testRTConstraints//Human_500nglysate_QTOF5600_tppSpstconsensus_lib_swathmsplit_out_pep_refResult_raw.txt";
		RTAligner align = createLocalAlignerFromResult(resultFile, 16, 17, 4, 6, 0, Double.MAX_VALUE);
		align.setPredictionMode(RTAligner.LOCAL);
		align.computeRegression();
		double[] interval = align.getRTDiffInterval();
		System.out.println("RT tolerance interval: " + interval[0] + " : " + interval[1]);
	}
	
	public static void main(String[] args){
		//testRTAligner();
		//testRunRTAlign();
		//testLocalRTAlign();
		testFilterMSPLITDIAwithRT();
	}
	
	
}
