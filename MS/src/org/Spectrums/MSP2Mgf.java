package org.Spectrums;
//wrapper for format conversion tools
public class MSP2Mgf {	
	public static void main(String[] args){
		String mgfFile = args[1];
		String mspFile = args[0];
		mgfFile = "..//mixture_linked//nistlib.mgf";
		mspFile = "../../Downloads/NIST_human_QTOF_2012-04-20.msp";
		if(args.length < 2){
			System.out.println("usage: java -Xmx1000M  -jar MSP2Mgf.jar <MSP file> <Mgf-out file>");	
		}else{
			System.out.println("Reading spectral library file");
			SpectrumLib lib = new SpectrumLib(mspFile, "MSP");
			System.out.println("Printing to MGF format");
			lib.printLibToFile(mgfFile, lib);
		}
	}
}
