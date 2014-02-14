package org.Spectrums;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Perform SVM classification of
 * @author Jian Wang
 *
 */
public class MixtureTDA extends TDAStat{
	public static int NO_MATCH = 0;
	public static int SINGLE_MATCH = 1;
	public static int MIX_MATCH = 2;
	private int protInd1=5;
	private int protInd2=5;
	private int scoreInd1 = 26;
	private int scoreInd2 = 26;
	private int rawScoreInd = 7;
	private double minScore = 0;
	private List<String> resultLines;
	private double threshold1;
	private double threshold2;
	private double tolerance = 10;
	
	public MixtureTDA(String resultFile){
		super(resultFile);
		this.resultLines = Utils.FileIOUtils.createListFromFile(this.resultFile);
	}
	public void filterByTDA(double FDR){
		filterByTDA(FDR, 10000);
	}
	
	public void filterByTDA(double FDR, double tolerance){
		List<Double> target1 = new ArrayList<Double>();
		List<Double> decoy1 = new ArrayList<Double>();
		for(Iterator<String> it = resultLines.iterator(); it.hasNext();){
			String current = it.next();
			String[] tokens = current.split("\\t");
			//System.out.println(current);
			this.tolerance = tolerance;
			if(checkResult(tokens)){
				//System.out.println(tokens[protInd1]);
				if(isDecoy(tokens[protInd1])){
					decoy1.add(Double.parseDouble(tokens[scoreInd1]));
				}else{
					target1.add(Double.parseDouble(tokens[scoreInd1]));
				}
			}	
		}
		this.threshold1 = getThreshold(target1, decoy1, FDR);
		
//		List<Double> target2 = new ArrayList<Double>();
//		List<Double> decoy2 = new ArrayList<Double>();
//		for(Iterator<String> it = resultLines.iterator(); it.hasNext();){
//			String current = it.next();
//			String[] tokens = current.split("\\t");
//			if(Double.parseDouble(tokens[rawScoreInd]) > this.minScore
//					&& (!isDecoy(tokens[protInd1]))
//					&& Double.parseDouble(tokens[scoreInd1]) > threshold1){
//				if(isDecoy(tokens[protInd2])){
//					decoy2.add(Double.parseDouble(tokens[scoreInd2]));
//				}else{
//					target2.add(Double.parseDouble(tokens[scoreInd2]));
//				}
//			}	
//		}
//		this.threshold2 = getThreshold(target2, decoy2, FDR);
		this.threshold2 = 100;
	}
	
	private boolean checkResult(String[] tokens){
		double mass1 = Double.parseDouble(tokens[2]);
		double mass2 = Double.parseDouble(tokens[6]);
		int charge = Integer.parseInt(tokens[3]);
		return (Double.parseDouble(tokens[rawScoreInd]) > this.minScore 
				&& !tokens[4].contains("Z")
				&& Mass.checkMassWithIso(mass1, mass2, tolerance, charge, 1, Mass.DIFF_PPM)
				&& tokens.length >= this.scoreInd2);
	}
	
	public void printOutFile(String outfile){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
			String header="#SpectrumFile\tScan#\tM/z\tCharge\tAnnotation\tProtein\tM/Z\tScore\tScore1\tScore2\tNormScore1\tNormScore2\t%Int\t%b_Pep\t%y_Pep\t%b_tag\t%y_tag\tbseries_pep\tyseries_pep\tbseries_tag\ty_series_tag\tmerror_pep\tmerrror_tag\n";
			int[] singleInds = new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
			int[] pairInds = new int[]{};
			int[] pairOutInds = new int[]{};
			out.write(header);
			for(Iterator<String> it = resultLines.iterator(); it.hasNext();){
				String current = it.next();
				String[] tokens = current.split("\\t");
				int match = getMatchClass(tokens);
				String[] outs = new String[singleInds.length];
				if(match != this.NO_MATCH){
					for(int i = 0; i < singleInds.length; i++){
						outs[i]=tokens[singleInds[i]-1];
					}
					outs[2]=outs[2].replaceAll("\\s+", ""); //it seems M-SPLIT print out an extra space at the end, remove it 
					
					if(match == this.MIX_MATCH){
						for(int i = 0; i < pairInds.length; i++){
							outs[pairOutInds[i]-1] = outs[pairOutInds[i]-1] 
								+ "!" + tokens[pairInds[i]-1]; 
						}
					}
				
					for(int i = 0; i < outs.length; i++){
						out.write(outs[i]+"\t");
					}
					out.write("\n");
				}
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
		
	private int getMatchClass(String[] tokens){
		int match = this.NO_MATCH;
		if(!checkResult(tokens)){
			return this.NO_MATCH;
		}
		if(Double.parseDouble(tokens[rawScoreInd]) > this.minScore
				&& (!isDecoy(tokens[protInd1]))
				&& Double.parseDouble(tokens[scoreInd1]) > threshold1){
			match = this.SINGLE_MATCH;
		}
		if(match == this.SINGLE_MATCH
				&& (!isDecoy(tokens[protInd2]))
				&& Double.parseDouble(tokens[scoreInd2]) > threshold2){
			match = this.MIX_MATCH;
		}
		return match;
	}
	
	
	public static void testTDAMixtre(){
		MixtureTDA mix = new MixtureTDA("../mixture_linked/MixDBv1.0/Human_heck_trypsin_mixdb_1_topHit_svmresult.txt");
		mix.filterByTDA(0.01);
	}
	
	
	public static void MixtureFilter(String resultFile, String outFile, double FDR, double tolerance, String svmtempDir){
		MixtureSVMClassify classify = new MixtureSVMClassify(resultFile, svmtempDir);
		classify.spectrumMatchClassify();
		String svmOut = classify.getSvmResultFile();
		MixtureTDA mixTDA = new MixtureTDA(svmOut);
		mixTDA.filterByTDA(FDR, tolerance);
		mixTDA.printOutFile(outFile);
		
	}
	
	public static void main(String[] args){
		if(args.length < 3 && args.length > 5){
			System.out.println("usage: java -Xmx1000M -jar MixDBFilter.jar <mixdb results> <output> <fdr> <precursor ppm_tolerance>");
			return;
		}
		//double fdr = 0.02;//Double.parseDouble(args[2]);
		//double tolerance = 25;//Double.parseDouble(args[3]);
		//String resultFile = "../mixture_linked/mxdb_sumo_outtop.txt";
		//String outFile = "../mixture_linked/mxdb_sumo_filteredout.txt";
	
		try{
			String resultFile = args[0];
			String outFile = args[1];
			double fdr = Double.parseDouble(args[2]);
			double tolerance = Double.parseDouble(args[3]);
				//double FDR = 0.01;
			String svmtempDir = System.getProperty("user.dir")+File.separator;
			if(args.length == 5){
				svmtempDir = args[4];
			}	
			MixtureFilter(resultFile, outFile, fdr, tolerance, svmtempDir);				
		}catch(Exception e){
				System.err.println(e.getMessage());
				e.printStackTrace();
				System.exit(-1);
		}
	}
}
