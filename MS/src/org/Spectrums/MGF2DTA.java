package org.Spectrums;
//convert mgf to dta file
//note in this format each dta file contain only one spectrum
public class MGF2DTA {
	public static void main(String[] args){
		String infile = args[0];
		String outdir = args[1];
		infile = "..\\mixture_linked\\human_heck_data\\data\\090121_NM_Trypsin_29.mgf";
		outdir = "..\\mixture_linked\\probidtree\\dtafiles\\heck_human_29";
		if(args.length < 2){
			System.out.println("usage: java -Xmx1000M  -jar MGF2DTA.jar <MGF file> <out dir>");	
		}else{
			System.out.println("Reading spectral library file");
			SpectrumLib lib = new SpectrumLib(infile, "MGF");
			lib.printStat(lib.NODETAIL);
			System.out.println("Printing to DTA format");
			lib.printLibToDTAFile(outdir, lib);
		}
	}
}
