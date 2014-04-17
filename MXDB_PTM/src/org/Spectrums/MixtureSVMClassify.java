package org.Spectrums;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Perform SVM classification for search result
 * Assume two svm model 
 * using svm_light
 * @author Jian Wang
 *
 */
public class MixtureSVMClassify {
	private static String currentBin= new File(MixtureSVMClassify.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent()
			+File.separator;
	private static String currentPath= System.getProperty("user.dir")+File.separator;
	//private String currentPath="../mixture_linked/MSPLIT_v1.0/";
	private String SVM_LIGHT_PATH= currentBin + "svm_light_windows" +
			"" + File.separator;	
	private String svmPath1 = SVM_LIGHT_PATH + "mxdb_PTM_sumo.model";
	private String svmPath2 = SVM_LIGHT_PATH + "mxdb_PTM_sumo.model";
	//private String svmPath2 = SVM_LIGHT_PATH + "msplit_stage2.model";
	private String resultFile;
	private String svmInFile;
	private String svmOutFile1;
	private String svmOutFile2;
	private String svmResultFile;
	private int beginFeatureInd=7;
	private int endFeatureInd=24;
	private int totalColumns = 25;
	private int[] skipIndices=new int[]{};
	
	public MixtureSVMClassify(String resultFile){
		this.resultFile = resultFile;
		initialize();
	}
	
	public MixtureSVMClassify(String resultFile, String classtempDir){
		this.resultFile = resultFile;
		this.currentPath = classtempDir;
		initialize();
	}
	
		

	public MixtureSVMClassify(String resultFile, String svmModel1, String svmModel2,
			int beginFeatureInd, int endFeatureInd, int totalColumns){
		this(resultFile, 
				MixtureSVMClassify.currentPath+File.separator,
				MixtureSVMClassify.currentBin + File.separator + "svm_light_windows" + File.separator,
				svmModel1, svmModel2, beginFeatureInd, endFeatureInd, totalColumns);
	}
	
	/**
	 * @param resultFile - raw result file
	 * @param classtempDir - optional a temporary dir to put temporary file for SVM, see constructor below
  	 * @param svmBinPath - optional path to svm binary 
	 * @param svmModel1  - name of first svm model file (assume in the binary path of svm)
	 * @param svmModel2  - name of second svm model file
	 * @param beginFeatureInd - begin index of result file to be used as feature
	 * @param endFeatureInd - end index of result file to be used as feature
	 * @param totalColumns
	 */
	public MixtureSVMClassify(String resultFile, String classtempDir, String svmBinPath, 
			String svmModel1, String svmModel2,
			int beginFeatureInd, int endFeatureInd, int totalColumns){
		this.resultFile = resultFile;
		this.currentPath = classtempDir;
		this.SVM_LIGHT_PATH = svmBinPath;
		this.svmPath1 = svmBinPath + svmModel1;
		this.svmPath2 = svmBinPath + svmModel2;
		this.beginFeatureInd = beginFeatureInd;
		this.endFeatureInd = endFeatureInd;
		this.totalColumns = totalColumns;
		initialize();
		//System.out.println("Current OS: " +  System.getProperty("os.name"));
		//System.out.println("path1: " + this.svmPath1);
	}
	
	
	//setup to run svm
	private void initialize(){
		this.svmResultFile = this.getSvmResultFile();
		this.svmInFile=currentPath + "temp_svmin.txt";
		this.svmOutFile1=currentPath + "temp_svmout1.txt";
		this.svmOutFile2=currentPath + "temp_svmout2.txt";
		this.svmResultFile=currentPath + "temp_svmresult.txt";
		System.out.println("out file: "  + this.svmResultFile);
	}
	
	/**
	 * Run the svm classification
	 * @return
	 */
	public String spectrumMatchClassify(){
		this.generateSVMInput();
		this.runSVMClassify();  
		this.getSVMResult();
		return this.svmResultFile;
	}
	
	private String getSVMResultFile(String inputFile){
		return inputFile + "_svmresult.txt";
	}
	
	private void generateSVMInput(){
		try{
			BufferedReader buff = new BufferedReader(new FileReader(this.resultFile));
			BufferedWriter out = new BufferedWriter(new FileWriter(this.svmInFile));
			String line = buff.readLine();
			//System.out.println("result file: " + this.resultFile);
			System.out.println("result file: " + this.svmInFile);
			while(line!=null){
				String[] tokens = line.split("\\t");
				//System.out.println("line: " + line);
				if(checkLine(line, tokens)){
					System.out.println("skipping " + line);
					line = buff.readLine();
					continue;
				}
				out.append("-1 ");
				int ind = 1;
				//System.out.println("skipIndice: " + skipIndices.length);
				for(int i = this.beginFeatureInd; i <= this.endFeatureInd; i++){
					boolean skip = false;
					for(int j = 0; j < this.skipIndices.length; j++){
						if(i == this.skipIndices[j]){
							skip=true;
						}
					}
					if(!skip){
						out.append((ind)+":"+tokens[i]+" ");
						ind++;
					}
				}
				out.append("\n");
				line = buff.readLine();
			}
			out.flush();
			out.close();
			buff.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	//check to see if the result line contains any error is comments that we should skipp
	private boolean checkLine(String line, String[] tokens){
		return line.contains("NaN") || tokens.length != this.totalColumns || line.startsWith("#");
	}
	
	private void runSVMClassify(){
		String cmd1 = this.SVM_LIGHT_PATH + "/svm_classify " + this.svmInFile + " " +  this.svmPath1 + " " + this.svmOutFile1;
		String cmd2 = this.SVM_LIGHT_PATH + "/svm_classify " + this.svmInFile + " " + this.svmPath2 + " " + this.svmOutFile2;
		//String[] cmd1 = new String[]{this.SVM_LIGHT_PATH + "/svm_classify.exe ", this.svmInFile + " " +  this.svmPath1, this.svmOutFile1};
		//String[] cmd2 = new String[]{this.SVM_LIGHT_PATH + "/svm_classify.exe ", this.svmInFile + " " +  this.svmPath2, this.svmOutFile1};
		try{
			System.out.println(cmd1);
			Process p1 = Runtime.getRuntime().exec(cmd1);
			gobbleStream(p1.getInputStream());
			gobbleStream(p1.getErrorStream());
			p1.waitFor();
			System.out.println(cmd2);
			Process p2 = Runtime.getRuntime().exec(cmd2);
			gobbleStream(p2.getInputStream());
			gobbleStream(p2.getErrorStream());
			p2.waitFor();
		}catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
	}
	
	//need to extract output from svm_classify otherwise process will not terminate
	private void gobbleStream(InputStream is){
		try{
			BufferedReader buff =new BufferedReader(new InputStreamReader(is));
			String line = buff.readLine();
			while(line != null){
				System.out.println(line);
				line = buff.readLine();
			}
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	private void getSVMResult(){
		try{
			BufferedReader buff = new BufferedReader(new FileReader(this.resultFile));
			BufferedReader svm1 = new BufferedReader(new FileReader(this.svmOutFile1));
			BufferedReader svm2 = new BufferedReader(new FileReader(this.svmOutFile2));
			BufferedWriter out = new BufferedWriter(new FileWriter(this.svmResultFile));
			String line = buff.readLine();
			while(line!=null){
				String[] tokens = line.split("\\t");
				if(checkLine(line, tokens)){
					line = buff.readLine();
					continue;
				}
				String svmvalue1 = svm1.readLine();
				String svmvalue2 = svm2.readLine();
				out.append(line + "\t" + svmvalue1 + "\t" + svmvalue2 + "\n");
				line = buff.readLine();
			}
			out.flush();
			out.close();
			buff.close();
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public String getCurrentPath() {
		return currentPath;
	}


	public void setCurrentPath(String currentPath) {
		this.currentPath = currentPath;
	}


	public String getSVM_LIGHT_PATH() {
		return SVM_LIGHT_PATH;
	}


	public void setSVM_LIGHT_PATH(String svm_light_path) {
		SVM_LIGHT_PATH = svm_light_path;
	}


	public String getSvmPath1() {
		return svmPath1;
	}


	public void setSvmPath1(String svmPath1) {
		this.svmPath1 = svmPath1;
	}


	public String getSvmPath2() {
		return svmPath2;
	}


	public void setSvmPath2(String svmPath2) {
		this.svmPath2 = svmPath2;
	}


	public String getResultFile() {
		return resultFile;
	}


	public void setResultFile(String resultFile) {
		this.resultFile = resultFile;
	}


	public String getSvmInFile() {
		return svmInFile;
	}


	public void setSvmInFile(String svmInFile) {
		this.svmInFile = svmInFile;
	}


	public String getSvmOutFile1() {
		return svmOutFile1;
	}


	public void setSvmOutFile1(String svmOutFile1) {
		this.svmOutFile1 = svmOutFile1;
	}


	public String getSvmOutFile2() {
		return svmOutFile2;
	}


	public void setSvmOutFile2(String svmOutFile2) {
		this.svmOutFile2 = svmOutFile2;
	}


	public String getSvmResultFile() {
		return svmResultFile;
	}


	public void setSvmResultFile(String svmResultFile) {
		this.svmResultFile = svmResultFile;
	}


	public int getEndFeatureInd() {
		return endFeatureInd;
	}


	public void setEndFeatureInd(int endFeatureInd) {
		this.endFeatureInd = endFeatureInd;
	}
	
	public int[] getSkipIndices() {
		return skipIndices;
	}

	public void setSkipIndices(int[] skipIndices) {
		this.skipIndices = skipIndices;
	}
	
	public int getTotalColumns() {
		return totalColumns;
	}

	
	public void setTotalColumns(int totalColumns) {
		this.totalColumns = totalColumns;
	}
		
	public static void testSVMClassify(){
		String resultFile = "../mixture_linked/t00";
		MixtureSVMClassify svmclassify = new MixtureSVMClassify(resultFile);
		svmclassify.spectrumMatchClassify();
	}
	
	public static void main(String[] args){
		testSVMClassify();
		System.out.println("finish running");
	}
}
