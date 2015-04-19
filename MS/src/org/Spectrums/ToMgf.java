package org.Spectrums;

///NOTE:   this is deprecated, should use the one in IO package
public class ToMgf {
	public static void MS22Mgf(String inFile, String outFile){
		System.out.println("Reading input library file");
		SpectrumLib lib = new SpectrumLib(inFile, "MS2");
		System.out.println("Printing to MGF format");
		lib.printLibToFile(outFile, lib);
	}
	
	public static void Sptxt2Mgf(String inFile, String outFile){
		System.out.println("Reading input library file");
		SpectrumLib lib = new SpectrumLib(inFile, "sptxt");
		System.out.println("Spectra list size: " + lib.getAllSpectrums().size());
		System.out.println("Printing to MGF format");
		lib.printLibToFile(outFile, lib);
	}
	
	public static void MSP2Mgf(String inFile, String outFile){
		System.out.println("Reading input library file");
		SpectrumLib lib = new SpectrumLib(inFile, "MSP");
		System.out.println("Spectra list size: " + lib.getAllSpectrums().size());
		System.out.println("Printing to MGF format");
		lib.printLibToFile(outFile, lib);
		
	}
	
	public static void main(String[] args){
		String inFile = args[0];
		String outFile = args[1];
		inFile = "..\\mixture_linked\\spectral_library\\NIST_mouse_consensus_final_true_lib_2014.sptxts";
		outFile = "..\\mixture_linked\\NIST_mouse_consensus_final_true_lib_2014.mgf";
		if(args.length < 2){
			System.out.println("usage: java -Xmx1000M  -jar ToMgf.jar <input file> <Mgf-out file>");	
		}else{
			Sptxt2Mgf(inFile, outFile);
			//MSP2Mgf(inFile, outFile);
		}
	}
}
