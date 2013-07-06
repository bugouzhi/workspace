package org.Spectrums;
//wrapper for format conversion tools
public class MSP2Mgf {	
	public static void main(String[] args){
		if(args.length < 2){
			System.out.println("usage: java -Xmx1000M  -jar MSP2Mgf.jar <MSP file> <Mgf-out file>");	
		}else{
			System.out.println("Reading spectral library file");
			SpectrumLib lib = new SpectrumLib(args[0], "MSP");
			System.out.println("Printing to MGF format");
			lib.printLibToFile(args[1], lib);
		}
	}
}
