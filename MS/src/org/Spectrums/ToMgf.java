package org.Spectrums;

public class ToMgf {
	public static void MS22Mgf(String inFile, String outFile){
		System.out.println("Reading input library file");
		SpectrumLib lib = new SpectrumLib(inFile, "MS2");
		System.out.println("Printing to MGF format");
		lib.printLibToFile(outFile, lib);
	}
	
	public static void Sptxt2Mgf(String inFile, String outFile){
		System.out.println("Reading input library file");
		SpectrumLib lib = new SpectrumLib(inFile, "splib");
		System.out.println("Printing to MGF format");
		lib.printLibToFile(outFile, lib);
	}
	
	public static void main(String[] args){
		String inFile = args[0];
		String outFile = args[1];
		inFile = "..\\mixture_linked\\spectral_library\\NIST_yeast_IT_v3.0_2009-05-04_7AA_decoy.sptxt";
		outFile = "..\\mixture_linked\\test.mgf";
		if(args.length < 2){
			System.out.println("usage: java -Xmx1000M  -jar ToMgf.jar <input file> <Mgf-out file>");	
		}else{
			Sptxt2Mgf(inFile, outFile);
		}
	}
}
