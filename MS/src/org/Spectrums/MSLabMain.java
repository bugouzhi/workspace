package org.Spectrums;
/**
 * This class contain usage info for using the various function in the whole 
 * MsLab package
 * @author Jian
 *
 */
public class MSLabMain {
	int status;
	public MSLabMain(){
		
	}
	
	public static void printUsage(){
		System.out.println("1) FeatureXMLParser <featurexml>");
		System.out.println("2) PeptideMassAnalysis <pep> <charge>");
		System.out.println("3) SpectrumUtil <spectrumDir> <annotationFile> <outfile> <7>");
	}
	
	public static void main(String[] args){
		printUsage();
	}
}
