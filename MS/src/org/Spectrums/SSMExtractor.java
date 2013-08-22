package org.Spectrums;

import java.util.List;

/**
 * Extract the matching peaks in a Spectrum-Spectrum Matches
 * from a spectral library search result and 
 * @author Jian Wang
 *
 */
public class SSMExtractor {
	private String querySpectrumFile;
	private String spectLibraryFile;
	private String searchResultFile;
	
	public SSMExtractor(String queryFile, String libFile, String results){
		this.querySpectrumFile = queryFile;
		this.spectLibraryFile = libFile;
		this.searchResultFile = results;
	}
	
	public void generateMatchedMGF(){
		List<String> results = Utils.FileIOUtils.createListFromFile(this.searchResultFile);
		MZXMLReader reader = new MZXMLReader(this.querySpectrumFile);
		SpectrumLib lib = new SpectrumLib(this.spectLibraryFile, "MGF");
		for(int i = 0; i < results.size(); i++){
			String line = results.get(i);
			String[] tokens = line.split("\\t");
			if(tokens.length < 6){
				continue;           //skipping header or other lines
			}
			Spectrum querySpect = reader.getSpectrum(Integer.parseInt(tokens[1]));
			//System.out.println(tokens[4] + "." + tokens[6]);
			//System.out.println(lib.getAllSpectrums().size());
			querySpect.windowFilterPeaks(15, 25);
			Spectrum libSpec = lib.getSpectra(tokens[4] + "." + Integer.parseInt(tokens[6])).get(0);
			Spectrum match = querySpect.project(libSpec, 0.05);
			match.spectrumName = querySpect.spectrumName;
			match.scanNumber = querySpect.scanNumber;
			match.parentMass = libSpec.parentMass;
			match.charge = libSpec.charge;
			match.peptide = libSpec.peptide;
			System.out.println(match);
			System.out.println();
		}
	}
	
	public static void main(String[] args){
		String query = "../mixture_linked/msdata/gringar/swath_development/1ugEcoli_400fmol_UPS2_swath_2012-12-11.mzXML";
		String lib = "../mixture_linked/ACG_swathdevelopment_P94_UPS_Ecoli_MSGFDB_IDs_uniqpeps_plusDecoy2_test.mgf";
		String result = "../mixture_linked/t00";
		SSMExtractor ssmE = new SSMExtractor(query, lib, result);
		ssmE.generateMatchedMGF();
	}
}
