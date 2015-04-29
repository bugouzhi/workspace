package Analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.Spectrums.ProteinIDExtractor;

/**
 * Compute the relative ratio of abundance from peptides from the same protein
 * assumption is that this ratio hold across data
 * @author Jian
 *
 */
public class PeptideRelativeRatio {
	Map<String, Double> pepRatioTable;
	int pepInd=0;
	int chargeInd=1;
	int[] abundanceInd= new int[]{1};//17, 19, 21, 23};
	boolean useCharge = false;
	
	public PeptideRelativeRatio(){
		
	}
	
	public void getPeptideRatioTable(String file){
		this.pepRatioTable = new HashMap<String, Double>();
		List<String> results = Utils.FileIOUtils.createListFromFile(file);
		//System.out.println("Number of lines: " + results.size());
		for(int i = 1; i < results.size(); i++){
			String[] tokens = results.get(i).split("\\t");
			double abundance = 0;
			int count=0;
			//System.out.println("line " + results.get(i));
			for(int j = 0; j < abundanceInd.length; j++){
				double value = Double.parseDouble(tokens[abundanceInd[j]]);
				if(value > 0){
					abundance+=value;
					count++;
				}
			}
			abundance /= count;
			//if(count > 3)
				//this.pepRatioTable.put(tokens[pepInd]+"@"+tokens[chargeInd], abundance);
				this.pepRatioTable.put(tokens[pepInd], abundance);
		}
		System.out.println("Got total peptides: " +  this.pepRatioTable.size());
	}
	
	public double getPredictedValue(String from, String to, double fromAbundance){
		if(this.pepRatioTable.containsKey(from) && this.pepRatioTable.containsKey(to)){
			double ratio = this.pepRatioTable.get(to) / this.pepRatioTable.get(from);
			//System.out.println("ratio is: " + ratio);
			return ratio*fromAbundance;
		}
		return -1.0;
	}
	
	public static void testPeptideRelativeRatio(){
		String mapFile = "..//mixture_linked//UPS_Human_quantstat.txt";
		PeptideRelativeRatio pepRatio = new PeptideRelativeRatio();
		pepRatio.getPeptideRatioTable(mapFile);
	}
	
	private static void getPepAbundanceMap(String file, Map<String, Set<String>> pepMap, Map<String, Double> abundanceMap){
		int pepInd= 1;
		int protInd = 0;
		int abundanceInd = 5;
		int chargeInd = 3;
		List<String> results = Utils.FileIOUtils.createListFromFile(file);
		
		for(int i = 1; i < results.size(); i++){ //skipping header
			String[] tokens = results.get(i).split("\\t");
			Set<String> peps = new HashSet<String>();
			//tokens[protInd] = "PROTEIN";
			if(pepMap.containsKey(tokens[protInd])){
				peps = pepMap.get(tokens[protInd]);
			}
			//peps.add(tokens[pepInd]+"@"+tokens[chargeInd]);
			peps.add(tokens[pepInd]);
			//pepMap.put(tokens[protInd], peps);
			pepMap.put(tokens[protInd], peps);
			//abundanceMap.put(tokens[pepInd]+"@"+tokens[chargeInd], Double.parseDouble(tokens[abundanceInd]));
			abundanceMap.put(tokens[pepInd], Double.parseDouble(tokens[abundanceInd]));
		}
	}
	
	public static void testPredictPepAbundance(){
		String mapFile = "..//mixture_linked//map.txt";
		String pepAbundanceFile = "..//mixture_linked//test.txt";
		String fasta = "..//mixture_linked//database//Human_uniprot-all.fasta";
		PeptideRelativeRatio pepRatio = new PeptideRelativeRatio();
		pepRatio.getPeptideRatioTable(mapFile);
		Map<String, Set<String>> pepMap = new HashMap<String, Set<String>>();
		Map<String, Double> abundanceMap = new HashMap<String, Double>();
		getPepAbundanceMap(pepAbundanceFile, pepMap, abundanceMap);
		ProteinIDExtractor protIDs = new ProteinIDExtractor(abundanceMap.keySet(), fasta);
		Map<String, List<String>> protMap = protIDs.getPeptideMap();
		for(Iterator<String> it = pepMap.keySet().iterator(); it.hasNext();){
			String prot = it.next();
			//Set<String> pep = pepMap.get(prot);
			List<String> pep = new ArrayList<String>();
			pep.addAll(pepMap.get(prot));
			for(int i = 0; i < pep.size(); i++){
				String pep1 = pep.get(i);
				List<Double> predicted = new ArrayList<Double>();
				for(int j = 0; j < pep.size(); j++){
					String pep2 = pep.get(j);
					if(i != j && protMap.get(pep1) != null && protMap.get(pep2) != null && (protMap.get(pep1).size() == 1 && protMap.get(pep2).size() == 1)){
						double predict21 = pepRatio.getPredictedValue(pep2, pep1, abundanceMap.get(pep2));
						if(predict21 > 0) predicted.add(predict21);
						double cv = Math.abs(predict21-abundanceMap.get(pep1)) / abundanceMap.get(pep1);
						if(predict21 > 0) System.out.println("Predicting abundance: "  + pep1+" & " + pep2 + "\t" + predict21 + "\t" + abundanceMap.get(pep1) + "\t" + abundanceMap.get(pep2) + "\tcv:\t" + cv);
					}
				}
				double predict1 = 0.0;
				for(int k = 0; k < predicted.size(); k++){
					predict1 += predicted.get(k);
				}
				if(predicted.size() > 0){
					Double[] arry = new Double[predicted.size()];
					predicted.toArray(arry);
					Arrays.sort(arry);
					predict1 = arry[(int)(arry.length/2.0)];
				}
				double cv = Math.abs(predict1-abundanceMap.get(pep1)) / abundanceMap.get(pep1);
				if(predict1 > 0) System.out.println("Consensus abundance: "  + pep1 + "\t" + predict1 + "\t" + abundanceMap.get(pep1) +  "\tcv:\t" + cv);	
			}
		}
		
	}
	
	public static void main(String[] args){
		//testPeptideRelativeRatio();
		testPredictPepAbundance();
	}
}
