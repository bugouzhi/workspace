

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Perform SVM classification for mixture spectrum
 * using svm_light
 * @author Jian Wang
 *
 */
public class MixtureSVMClassify {
	private String currentPath=System.getProperty("user.dir")+File.separator;
	//private String currentPath="../mixture_linked/MSPLIT_v1.0/";
	private String SVM_LIGHT_PATH= getDefaultSVMPath();
	//private String svmPath1 = SVM_LIGHT_PATH + "mixdb_stage1.model";
	//private String svmPath1 = SVM_LIGHT_PATH + "human_simmixmixdb_stage1.model";
	private String svmPath1 = SVM_LIGHT_PATH + "msplit_stage1.model";
	//private String svmPath2 = SVM_LIGHT_PATH + "mixdb_stage2.model";
	//private String svmPath2 = SVM_LIGHT_PATH + "human_simmixmixdb_stage2.model";
	private String svmPath2 = SVM_LIGHT_PATH + "msplit_stage2.model";
	private String resultFile=currentPath;
	private String svmInFile=currentPath + "temp_svmin.txt";
	private String svmOutFile1=currentPath + "temp_svmout1.txt";
	private String svmOutFile2=currentPath + "temp_svmout2.txt";
	private String svmResultFile=currentPath + "temp_svmresult.txt";
	private int beginFeatureInd=8;//10;
	private int endFeatureInd=22;//25;
	private int rawScoreInd = 8;//10;
	private double minScore = 0.4;//32;
	private int totalColumns = 29;//28;
	
	public static String getDefaultSVMPath(){
		String bin = new File(MixtureSVMClassify.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent()
				+File.separator;
		return bin + "svm_light_linux" + "" + File.separator;
	}
	
	public MixtureSVMClassify(String resultFile){
		this.resultFile = resultFile;
		this.svmResultFile = this.getSVMResultFile(this.resultFile);
		//System.out.println("out file: "  + this.svmResultFile);
		//this.currentPath = System.getProperty("user.dir");
		//System.out.println("currentpath is: " + currentPath);
		//this.currentPath = "C:\\Documents and Settings\\Jian Wang\\workspace\\mixture_linked/MSPLIT";
		//this.svmInFile = currentPath +"/temp_svmin.txt";
		//this.svmOutFile1 = currentPath + "/temp_svmout1.txt";
		//this.svmOutFile2 = currentPath + "/temp_svmout2.txt";
		//this.svmResultFile = currentPath + "/temp_svmresult.txt";
	}
	
	public MixtureSVMClassify(String resultFile, String classtempDir){
		this.resultFile = resultFile;
		this.svmResultFile = this.getSVMResultFile(this.resultFile);
		//System.out.println("out file: "  + this.svmResultFile);
		this.currentPath = classtempDir;
		System.out.println("currentpath is: " + this.currentPath);
		//this.currentPath = "C:\\Documents and Settings\\Jian Wang\\workspace\\mixture_linked/MSPLIT";
		this.svmInFile = currentPath +"/temp_svmin.txt";
		this.svmOutFile1 = currentPath + "/temp_svmout1.txt";
		this.svmOutFile2 = currentPath + "/temp_svmout2.txt";
		this.svmResultFile = currentPath + "/temp_svmresult.txt";
	}

	
	public MixtureSVMClassify(String resultFile, String classtempDir, String svmPath){
		this.resultFile = resultFile;
		this.SVM_LIGHT_PATH = svmPath;
		this.svmResultFile = this.getSVMResultFile(this.resultFile);
		//System.out.println("out file: "  + this.svmResultFile);
		this.currentPath = classtempDir;
		System.out.println("currentpath is: " + this.currentPath);
		//this.currentPath = "C:\\Documents and Settings\\Jian Wang\\workspace\\mixture_linked/MSPLIT";
		this.svmInFile = currentPath +"/temp_svmin.txt";
		this.svmOutFile1 = currentPath + "/temp_svmout1.txt";
		this.svmOutFile2 = currentPath + "/temp_svmout2.txt";
		this.svmResultFile = currentPath + "/temp_svmresult.txt";
	}
	
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
				if(line.contains("NaN") || tokens.length != this.totalColumns || line.contains("#") 
						|| Double.parseDouble(tokens[this.rawScoreInd]) < this.minScore){
					line = buff.readLine();
					continue;
				}
				out.append("-1 ");
				for(int i = this.beginFeatureInd; i <= this.endFeatureInd; i++){
					out.append((i-this.beginFeatureInd+1)+":"+tokens[i]+" ");
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
	
	private void runSVMClassify(){
		String cmd1 = this.SVM_LIGHT_PATH + "/svm_classify " + this.svmInFile + " " +  this.svmPath1 + " " + this.svmOutFile1;
		String cmd2 = this.SVM_LIGHT_PATH + "/svm_classify " + this.svmInFile + " " + this.svmPath2 + " " + this.svmOutFile2;
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
				if(line.contains("NaN") || tokens.length != this.totalColumns|| line.startsWith("#")
						|| Double.parseDouble(tokens[this.rawScoreInd]) < this.minScore){
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


	public int getRawScoreInd() {
		return rawScoreInd;
	}


	public void setRawScoreInd(int rawScoreInd) {
		this.rawScoreInd = rawScoreInd;
	}


	public double getMinScore() {
		return minScore;
	}


	public void setMinScore(double minScore) {
		this.minScore = minScore;
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
