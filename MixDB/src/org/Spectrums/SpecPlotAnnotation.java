package org.Spectrums;

import java.util.List;

import Utils.FileIOUtils;

/**
 * Generate annotation file for specplot
 * @author Jian Wang
 *
 */
public class SpecPlotAnnotation {
	private static MZXMLReader reader = null;
	public static boolean DEBUG = false;
	public static void getLinkedSpectPlot(String spectrumLibFile, String annotationFile, double tolerance, String outdir){
		List<String> lines = FileIOUtils.createListFromFile(annotationFile);
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		for(int i = 0; i < lines.size(); i++){
			String[] tokens = lines.get(i).split("\\s+");
			int scan = Integer.parseInt(tokens[2]);
			String peptide = tokens[3];
			int position1 = Integer.parseInt(tokens[4]);
			int position2 = Integer.parseInt(tokens[5]);
			double linkerMass = Double.parseDouble(tokens[6]);
			int charge = -1;
			if(tokens.length == 8){
				charge = Integer.parseInt(tokens[7]);
			}
			getLinkedSpectPlot(spectrumLibFile, scan, peptide, position1, position2, linkerMass, 0.5, charge, outdir);
		}
	}
	
	public static void getLinkedSpectPlot(String spectrumLibFile, int scan, String peptide, int position1, int position2, double linkermass, double tolerance, int charge, String outdir){
		if(reader == null){
			reader = new MZXMLReader(spectrumLibFile);
		}
		Spectrum s = reader.getSpectrum(scan);
		//s.windowFilterPeaks(1, 1.0);
		s.windowFilterPeaks(8, 25);
		//int charge = Integer.parseInt(tokens[5]);
		if(charge < 0){
			charge = s.charge;
		}
		if(charge < 0 ){
			System.out.println("warnin: no charge information");
			return;
		}
		TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		Mass.DSSLINKER_MASS = linkermass;
		System.out.println("peptide is: " + peptide + " charge is : " + s.charge);
		LinkedPeptide lp = new LinkedPeptide(peptide, charge, position1, position2);
		System.out.println("peptide is: " + lp + "\t" + lp.getCharge() + "\t" + lp.getParentmass());
		String outfile = outdir+"spectPlot_"+scan+"_"+peptide+".txt";
		TheoreticalSpectrum linkedSpect = new TheoreticalSpectrum(lp, lp.getCharge(), false);
		s.computePeakRank();
		linkedSpect.analyzeMixtureAnnotation(s, lp.peptides[0].getPeptide(), lp.peptides[1].getPeptide(), tolerance);
		linkedSpect.printSpectPlotAnnotationFile(s, outfile, tolerance);
	}
	
	
	public static void main(String[] args){
		args = new String[4];
		args[0] = "..\\mixture_linked\\linked_peptide_library\\UWash/BSA_DSS/030909_SCX_300mM.mzXML";
		args[1] = "..\\mixture_linked\\testAnnotation.txt";
		args[2] = "0.5";
		args[3] = "..\\mixture_linked\\specplotAnnotation\\";
		if(args.length == 9){
			getLinkedSpectPlot(args[0], Integer.parseInt(args[1]), args[2], Integer.parseInt(args[3]), 
					Integer.parseInt(args[4]), Double.parseDouble(args[5]), Double.parseDouble(args[6]), Integer.parseInt(args[7]), args[8]);
		}else if(args.length == 4){
			getLinkedSpectPlot(args[0], args[1], Double.parseDouble(args[2]), args[3]);
		}else{
			System.out.println("Usage: java SpecPlotAnnotation <spectrumLibFile>  <scan> <peptide> <position1> <position2> <linkermass> <tolerance> < charge> <outdir>");
			System.out.println("Usage: java SpecPlotAnnotation <spectrumLibFile> <annotationFile> <tolerance> <outdir>");
		}
	}

}
