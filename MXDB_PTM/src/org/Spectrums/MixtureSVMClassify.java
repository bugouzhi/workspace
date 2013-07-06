package org.Spectrums;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
	private static String SVM_LIGHT_PATH="../mixture_linked/MSPLIT/svm_light_windows";
	private String resultFile;
	private String svmInFile="../mixture_linked/temp_svmin.txt";
	private String svmOutFile1="../mixture_linked/temp_svmout1.txt";
	private String svmOutFile2="../mixture_linked/temp_svmout2.txt";
	private String svmResultFile="../mixture_linked/temp_svmresult.txt";
	private int beginFeatureInd=8;
	private int endFeatureInd=21;
	
	public MixtureSVMClassify(String resultFile){
		this.resultFile = resultFile;
	}
	
	public void spectrumMatchClassify(){
		this.generateSVMInput();
		this.runSVMClassify();
		this.getSVMResult();
	}
	
	private void generateSVMInput(){
		try{
			BufferedReader buff = new BufferedReader(new FileReader(this.resultFile));
			BufferedWriter out = new BufferedWriter(new FileWriter(this.svmInFile));
			String line = buff.readLine();
			//System.out.println("result file: " + this.resultFile);
			while(line!=null){
				String[] tokens = line.split("\\t");
				//System.out.println("line: " + line);
				if(line.contains("NaN") || tokens.length != 26 || line.contains("#")){
					line = buff.readLine();
					continue;
				}
				out.append("1 ");
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
		String cmd1 = this.SVM_LIGHT_PATH + "/svm_classify.exe " + this.svmInFile + " " + this.SVM_LIGHT_PATH + "/msplit_stage1.model " + this.svmOutFile1;
		String cmd2 = this.SVM_LIGHT_PATH + "/svm_classify.exe " + this.svmInFile + " " + this.SVM_LIGHT_PATH + "/msplit_stage2.model " + this.svmOutFile2;
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
				if(line.contains("NaN") || tokens.length != 26 || line.contains("#")){
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
		
	public static void testSVMClassify(){
		String resultFile = "../mixture_linked/MSPLIT/test.txt";
		MixtureSVMClassify svmclassify = new MixtureSVMClassify(resultFile);
		svmclassify.spectrumMatchClassify();
	}
	
	public static void main(String[] args){
		testSVMClassify();
		System.out.println("finish running");
	}
}
