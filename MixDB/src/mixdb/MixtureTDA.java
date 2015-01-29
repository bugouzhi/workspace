package mixdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Perform SVM classification of
 * @author Jian Wang
 *
 */
public class MixtureTDA {
	public static int NO_MATCH = 0;
	public static int SINGLE_MATCH = 1;
	public static int MIX_MATCH = 2;
	private String resultFile;
	private int protInd1=7;
	private int protInd2=8;
	private int scoreInd1 = 27;
	private int scoreInd2 = 28;
	private int rawScoreInd = 9;
	private double minScore = 30;
	private List<String> results;
	private double threshold1;
	private double threshold2;
	private double sortMode = 1;
	
	public MixtureTDA(String resultFile){
		this.resultFile = resultFile;
		this.results = Utils.FileIOUtils.createListFromFile(this.resultFile);
	}
	
	public void filterByTDA(double FDR){
		List<Double> target1 = new ArrayList<Double>();
		List<Double> decoy1 = new ArrayList<Double>();
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String current = it.next();
			String[] tokens = current.split("\\t");
			//System.out.println("current line is : " + current + "\t" + tokens.length);
			if(tokens.length < scoreInd2){
				continue;
			}
			if(Double.parseDouble(tokens[rawScoreInd]) > this.minScore){
				//System.out.println("Score is: " + tokens[scoreInd1]);
				if(isDecoy(tokens[protInd1])){
					decoy1.add(this.sortMode*Double.parseDouble(tokens[scoreInd1]));
				}else{
					target1.add(this.sortMode*Double.parseDouble(tokens[scoreInd1]));
				}
			}	
		}
		this.threshold1 = getThreshold(target1, decoy1, FDR);
		
		List<Double> target2 = new ArrayList<Double>();
		List<Double> decoy2 = new ArrayList<Double>();
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String current = it.next();
			String[] tokens = current.split("\\t");
			if(tokens.length < scoreInd2){
				continue;
			}
			if(Double.parseDouble(tokens[rawScoreInd]) > this.minScore
					&& (!isDecoy(tokens[protInd1]))
					&& this.sortMode*Double.parseDouble(tokens[scoreInd1]) > threshold1){
				if(isDecoy(tokens[protInd2])){
					decoy2.add(this.sortMode*Double.parseDouble(tokens[scoreInd2]));
				}else{
					target2.add(this.sortMode*Double.parseDouble(tokens[scoreInd2]));
				}
			}	
		}
		this.threshold2 = getThreshold(target2, decoy2, FDR);
	}
	
	public void printOutFile(String outfile){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
			String header="#SpectrumFile\tScan#\tAnnotation\tProtein\tscore\tscore-1pep\texpl-Int\tb-fract\ty-fract\tb-series\ty-series\tsvm1-score\tsvm2-score\n";
			int[] singleInds = new int[]{1,2,6,8,10,11,15,16,17,20,21,28,29};
			int[] pairInds = new int[]{7,9,12,18,19,22,23};
			int[] pairOutInds = new int[]{3,4,6,8,9,10,11};
			out.write(header);
			for(Iterator<String> it = results.iterator(); it.hasNext();){
				String current = it.next();
				String[] tokens = current.split("\\t");
				tokens[5]=tokens[5].replaceAll("[0-9.]", "");
				tokens[6]=tokens[6].replaceAll("[0-9.]", "");
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
	
	
	public void printOutFile(String outfile, String header, int[] singleInds, int[] pairInds, int[] pairOutInds){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
			out.write(header);
			for(Iterator<String> it = results.iterator(); it.hasNext();){
				String current = it.next();
				String[] tokens = current.split("\\t");
				if(tokens.length < this.scoreInd2){
					continue;
				}
				tokens[5]=tokens[5].replaceAll("[0-9.]", "");
				tokens[6]=tokens[6].replaceAll("[0-9.]", "");
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
		if(Double.parseDouble(tokens[rawScoreInd]) >= this.minScore
				&& (!isDecoy(tokens[protInd1]))
				&& this.sortMode*Double.parseDouble(tokens[scoreInd1]) > threshold1){
			match = this.SINGLE_MATCH;
		}
		if(match == this.SINGLE_MATCH
				&& (!isDecoy(tokens[protInd2]))
				&& this.sortMode*Double.parseDouble(tokens[scoreInd2]) >= threshold2){
			match = this.MIX_MATCH;
		}
		return match;
	}
	
	
	private double getThreshold(List<Double> target, List<Double> decoy, double fdr){
		int totalTarget = target.size();
		int totalDecoy = decoy.size();
		System.out.println("Total numbers of targets: " + target.size());
		System.out.println("Total numbers of decoys: " + decoy.size());
		Collections.sort(target);
		Collections.sort(decoy);
		int i = 0, j = 0; 
		while(i < target.size() && j < decoy.size()){
			//System.out.println("target: " + (target.size()-i) +"\t" + "decoy: " + (decoy.size()-j) + "\t" + target.get(i) + "\t" + decoy.get(j));
			if(target.get(i) < decoy.get(j)){
				if((double)(totalDecoy - j ) / (double)(totalTarget - i ) < fdr){
					break;
				}
				i++;
			}else{
				j++;
			}
		}
		
		if(i == target.size() && i > 0){
			i=i-1;
		}
		
		System.out.println("threshold is: " + target.get(i) + " Accepted targets: " + (totalTarget - i));
		return target.get(i);
	}
	
	
	private boolean isDecoy(String prot){
		return prot.contains("DECOY_") || prot.startsWith("X_");
	}
	
	public static void testTDAMixtre(){
		MixtureTDA mix = new MixtureTDA("../mixture_linked/MixDBv1.0/Human_heck_trypsin_mixdb_1_topHit_svmresult.txt");
		mix.filterByTDA(0.01);
	}
	public static void MixtureFilter(String resultFile, String outFile, double FDR, String tempDir){
		MixtureSVMClassify classify = new MixtureSVMClassify(resultFile, tempDir);
		classify.spectrumMatchClassify();
		String svmOut = classify.getSvmResultFile();
		MixtureTDA mixTDA = new MixtureTDA(svmOut);
		mixTDA.filterByTDA(FDR);
		mixTDA.printOutFile(outFile);
	}
	
	public static void MixGFfilter(String resultFile, String outFile, double FDR){
		String header="#SpectrumFile\tScan#\tAnnotation\tProtein\tscore\tscore-1pep\texpl-Int\tb-fract\ty-fract\tb-series\ty-series\tSingle-Prob\tCond-Prob\n";
		int[] singleInds = new int[]{1,2,6,8,10,11,15,16,17,20,21,32,33};
		int[] pairInds = new int[]{7,9,12,18,19,22,23};
		int[] pairOutInds = new int[]{3,4,6,8,9,10,11};
		MixtureTDA mixTDA = new MixtureTDA(resultFile);
		mixTDA.minScore = 0;
		mixTDA.scoreInd1 = 32;
		mixTDA.scoreInd2 = 33;
		mixTDA.sortMode = -1;
		mixTDA.filterByTDA(FDR);
		mixTDA.printOutFile(outFile, header, singleInds, pairInds, pairOutInds);
	}
	
	public static void parseInputs(String inputFile){
		Map<String, String> table = Utils.FileIOUtils.createTableFromFile(inputFile, 0, 1);
	}
	
	public static void main(String[] args){
		if(args.length != 3 && args.length != 4){
			System.out.println("usage: java -Xmx1000M -jar MixDB/GFFilter.jar <mixdb results> <output> <fdr>");
			return;
		}
		//String resultFile = "../mixture_linked/2MixturesApha1.0_mixgfout.txt";
		//String outFile = "../mixture_linked/MixtureFiltered.test.txt";
		try{
			String resultFile = args[0];
			String outFile = args[1];
			double FDR = Double.parseDouble(args[2]);
			//double FDR = 0.02;
			String svmtempDir = System.getProperty("user.dir")+File.separator;
			if(args.length == 4){
				svmtempDir = args[3];
			}
			//MixtureFilter(resultFile, outFile, FDR, svmtempDir);
			MixGFfilter(resultFile, outFile, FDR);
		}catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
