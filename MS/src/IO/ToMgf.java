package IO;

import org.Spectrums.SpectrumLib;
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
		System.out.println("Spectra list size: " + lib.getAllSpectrums().size());
		System.out.println("Printing to MGF format");
		lib.printLibToFile(outFile, lib);
	}
	
	public static void ToMgf(String inFile, String outFile, String inFormat){
		System.out.println("Reading input library file");
		SpectrumLib lib = new SpectrumLib(inFile, inFormat);
		//lib.removeModSpectra();
		System.out.println("Spectra list size: " + lib.getAllSpectrums().size());
		System.out.println("Printing to MGF format");
		lib.printLibToFile(outFile, lib);
	}
	
	public static void main(String[] args){
		String inFile = args[0];
		String outFile = args[1];
		inFile = "..\\mixture_linked\\HSA_NIST_speclib\\hsa20140422.sptxt";
		outFile = "..\\mixture_linked\\HSA_NIST_speclib\\test.mgf";
		if(args.length < 3){
			System.out.println("usage: java -Xmx1000M  -jar ToMgf.jar <input file> <Mgf-out file> <input format>");	
		}else{
			ToMgf(inFile, outFile, "sptxt");
		}
	}
}
