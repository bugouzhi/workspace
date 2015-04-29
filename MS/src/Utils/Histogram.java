package Utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import UI.CommandLineParser;

/**
 * A simple class to generate text base histogram
 *  * @author Jian
 *
 */
public class Histogram {
	private List<int[]> histCounts;
	double[][] boundaries;
	private int[] indexes;
	private String resultFile;
	public Histogram(String resultFile){
		this.resultFile = resultFile;
		this.histCounts = new ArrayList<int[]>();
	}
	
	/**
	 * This method stream data from file, can scale to rather big data
	 * @param ind
	 * @param boundaries
	 * @return
	 */
	public List<int[]> generateHistFromFile(int[] ind, double[][] boundaries){
		BufferedReader reader = FileIOUtils.createReaderFromFile(this.resultFile);
		this.histCounts.clear();
		this.boundaries = boundaries;
		if(ind.length != boundaries.length){
			throw new IllegalArgumentException("Index and boundaries do no have the same size");
		}
		for(int i = 0; i < ind.length; i++){
			histCounts.add(new int[boundaries[i].length-1]);
		}
		try{
			String line = reader.readLine();
			while(line != null){
				if(line.startsWith("#")){
					line = reader.readLine();
					continue;
				}
				String[] tokens = line.split("\\t");
				for(int i = 0; i < ind.length; i++){
					double value = Double.parseDouble(tokens[ind[i]]);
					//System.out.println("data is " + value + "\tind is: " + ind[i]);
					int[] histCount = histCounts.get(i);
					histCount[ArrayUtils.getIntervalIndex(value, boundaries[i])]++;
				}
				//System.out.println("line is " + line);
				line = reader.readLine();
			}
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		return getHistCounts();
	}
	
	public List<int[]> generateHistFromData(double[][] data, double[][] boundaries){
		for(int i = 0; i < data.length; i++){
			int[] histCount = ArrayUtils.hist1D(data[i], boundaries[i]);
			this.histCounts.add(histCount);
		}
		return getHistCounts();
	}
	
	
	
	public List<int[]> getHistCounts() {
		List<int[]> retrn = new ArrayList<int[]>();
		retrn.addAll(this.histCounts);
		return retrn;
	}

	public void setHistCounts(List<int[]> histCounts) {
		this.histCounts = histCounts;
	}

	public void displayHist(String outFile){
		displayHist(outFile, null);
	}
	public void displayHist(String outFile, List<String> titles){
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		for(int i = 0; i < histCounts.size(); i++){
			if(titles != null){
				out.println(titles.get(i));
			}
			out.print(ArrayUtils.displayHist1D(histCounts.get(i), boundaries[i]));
			out.println();
		}
	}
	
	public static void testMultiHistFromFile(String resultFile, String paramFile, String outFile){
		Histogram hist = new Histogram(resultFile);
		List<String> histParams = Utils.FileIOUtils.createListFromFile(paramFile);
		int[] inds = new int[histParams.size()];
		double[][] boundaries = new double[histParams.size()][];
		List<String> titles = new ArrayList(histParams.size());
		for(int i = 0; i < histParams.size(); i++){
			String[] tokens = histParams.get(i).split(":");
			if(tokens.length < 3){
				throw new IllegalArgumentException("Histogram parameters format not valid, should be:  <Title>:<index>:<boundary1>,<boundary2>,.....");
			}
			titles.add(tokens[0]);
			inds[i] = Integer.parseInt(tokens[1]);
			boundaries[i] = StringUtils.toDoubleArray(tokens[2].split(","));
		}
		hist.generateHistFromFile(inds, boundaries);
		hist.displayHist(outFile, titles);
	}
	
	
	
	public static void testHistFromFile(String args[]){
		String resultFile = "..//mixture_linked//createlib.log";
		String index = "21";
		String boundaries = "0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,"+Integer.MAX_VALUE;
		String out = "";
		runHistFromFile(new String[]{resultFile, index, boundaries, out});
	}
	
	public static void runHistFromFile(String args[]){
		CommandLineParser parser = new CommandLineParser(args);
		String resultFile = parser.getString(0);
		//System.out.println("result file " + resultFile);
		if(args.length == 3){
			testMultiHistFromFile(resultFile, parser.getString(1), parser.getString(2));
		}else{
			int index = parser.getInteger(1);
			double[] boundaries = parser.getDoubles(2);
			Histogram hist = new Histogram(resultFile);
			hist.generateHistFromFile(new int[]{index}, new double[][]{boundaries});
			if(parser.isAvailable(3)){
				hist.displayHist(parser.getString(3));
			}else{
				hist.displayHist("");
			}
		}
	}
	
	public static void main(String[] args){
		//testHistFromFile(args);
		runHistFromFile(args);
	}
	
	
	
}
