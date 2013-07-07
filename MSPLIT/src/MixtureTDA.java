

import java.io.BufferedWriter;
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
public class MixtureTDA {
	public static int NO_MATCH = 0;
	public static int SINGLE_MATCH = 1;
	public static int MIX_MATCH = 2;
	private String resultFile;
	private int protInd1=4;//7;
	private int protInd2=5;//8;
	private int scoreInd1 = 27;
	private int scoreInd2 = 28;
	private int rawScoreInd = 8;//9;
	private double minScore = 0.4;//0;
	private List<String> results;
	private double threshold1;
	private double threshold2;
	
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
			if(Double.parseDouble(tokens[rawScoreInd]) > this.minScore){
				if(isDecoy(tokens[protInd1])){
					decoy1.add(Double.parseDouble(tokens[scoreInd1]));
				}else{
					target1.add(Double.parseDouble(tokens[scoreInd1]));
				}
			}	
		}
		this.threshold1 = getThreshold(target1, decoy1, FDR);
		
		List<Double> target2 = new ArrayList<Double>();
		List<Double> decoy2 = new ArrayList<Double>();
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String current = it.next();
			String[] tokens = current.split("\\t");
			if(Double.parseDouble(tokens[rawScoreInd]) > this.minScore
					&& (!isDecoy(tokens[protInd1]))
					&& Double.parseDouble(tokens[scoreInd1]) > threshold1){
				if(isDecoy(tokens[protInd2])){
					decoy2.add(Double.parseDouble(tokens[scoreInd2]));
				}else{
					target2.add(Double.parseDouble(tokens[scoreInd2]));
				}
			}	
		}
		this.threshold2 = getThreshold(target2, decoy2, FDR);
	}
	
	public void printOutFile(String outfile){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
			String header="#SpectrumFile\tScan#\tAnnotation\tProtein\tCharge\tcosine(M, A+B)\tcosine(M,A)\tcosine(A,B)\talpha\tres-alpha\t#Peak-0.85Intensity\tsimBias(M,A+B)\tsimBias(A)\tprojCos(M,A+B)\tprojCos(M,A)\tmeanCos\tmeanDeltaCos\tPrecusor(M)\tPrecursor(A)\tspectrumIndex\tsvm1-score\tsvm2-score\n";
			int[] singleInds = new int[]{1,2,3,5,7,9,10,12,13,14,15,16,17,19,21,22,23,24,25,27,28,29};
			int[] pairInds = new int[]{4,6,8,11,18,21,26};
			int[] pairOutInds = new int[]{3,4,5,7,13,14,19};
			out.write(header);
			for(Iterator<String> it = results.iterator(); it.hasNext();){
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
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
		
	private int getMatchClass(String[] tokens){
		int match = this.NO_MATCH;
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
			
			if(target.get(i) <= decoy.get(j)){
				i++;
			}else{
				j++;
			}
			if((double)(totalDecoy - j ) / (double)(totalTarget - i ) < fdr){
				break;
			}
		}
		System.out.println("threshold is: " + target.get(i-1) + " Accepted targets: " + (totalTarget - i));
		return target.get(i);
	}
	
	
	private boolean isDecoy(String prot){
		return prot.contains("DECOY_") || prot.startsWith("X_");
	}
	
	public static void testTDAMixtre(){
		MixtureTDA mix = new MixtureTDA("../mixture_linked/MixDBv1.0/Human_heck_trypsin_mixdb_1_topHit_svmresult.txt");
		mix.filterByTDA(0.01);
	}
	public static void MixtureFilter(String resultFile, String outFile, double FDR){
		MixtureSVMClassify classify = new MixtureSVMClassify(resultFile);
		classify.spectrumMatchClassify();
		String svmOut = classify.getSvmResultFile();
		MixtureTDA mixTDA = new MixtureTDA(svmOut);
		mixTDA.filterByTDA(FDR);
		mixTDA.printOutFile(resultFile +".filtered.txt");
		
	}
	
	public static void main(String[] args){
		String resultFile = args[0];
		String outFile = args[1];
		double fDR = Double.parseDouble(args[2]);
		//resultFile = "../mixture_linked/MSPLIT_v1.0/40fmol_UPS2_1ugEcoli_MSPLIT_testsearch.txt";
		//fDR = 0.01;
		MixtureFilter(resultFile, outFile, fDR);
	}
}
