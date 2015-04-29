package Misc;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.Spectrums.ProteinIDExtractor;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumLib;



/**
 * Various scripts use for processing spectral library
 * @author Jian
 *
 */
public class SpectrumLibScripts {
	/**
	 * Map protein entry to SNP spectral lib
	 * The main purpose is to assign proper protein name so peptide with SNP will get
	 * an SNP entry from the database and peptides with no SNP will map to protein with 
	 * regular entries (when more than one entry available, map to the first one)
	 * @param specLibFile
	 * @param fastaFile
	 */
	public static void getSNPProteinEntries(String specLibFile, String fastaFile){
		SpectrumLib.USECHARGE_DEFAULT = true;
		SpectrumLib lib = new SpectrumLib(specLibFile, "MGF");
		Set<String> peps = new HashSet(lib.getSpectrumLibrary().keySet());
		peps.addAll(lib.getSpectrumLibrary().keySet());
		ProteinIDExtractor protIDs = new ProteinIDExtractor(peps, fastaFile);
		PrintStream out = Utils.FileIOUtils.getOutStream(Utils.FileIOUtils.stripExtension(specLibFile)+"_withSNPProt.mgf");
		for(int i = 0; i < lib.getAllSpectrums().size(); i++){
			Spectrum s = (Spectrum) lib.getAllSpectrums().get(i);
			if(s.spectrumName.contains("DECOY")){
				continue;
			}
			List<String> prots = protIDs.getProteins(s.peptide);
			if(prots.size() > 0){
				s.protein = prots.get(0);
			}
			for(int j = 0; j < prots.size(); j++){
				if(!isSNP(prots.get(j))){
					s.protein = prots.get(j);
					break;
				}
			}
			out.println(s);
		}
	}
	
	private static boolean isSNP(String prot){
		return prot.contains("dbSNP");
	}
	
	public static void main(String[] args){
		String lib = "../mixture_linked/SunJun_Privacy_Data/Human_SNP_lib.mgf";
		String fastaFile = "../mixture_linked/SunJun_Privacy_Data/160_Glyco_MSGF_FDR1.withdbSNP.pep.fasta";
		getSNPProteinEntries(lib, fastaFile);
		
	}
}
