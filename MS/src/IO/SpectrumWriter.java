package IO;

import java.io.PrintStream;
import java.util.List;

import org.Spectrums.Peak;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumLib;

/**
 * Write out MS spectrum in specified format
 * @author Jian
 *
 */
public class SpectrumWriter {
	private String outFormat;
	public void printSpectraToFile(List<Spectrum> specList, String outFile){
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			out.println(this.toMSPString(s));
		}
	}
	
	public String toMSPString(Spectrum s){
		String RT = "0";
		if(s.spectrumName.split("\\s+").length > 5 && s.spectrumName.contains("Retention Time")){  //check if retention time is in lib
		//	System.out.println("SpectruMName: " + cand.spectrumName + "\t" + cand.spectrumName.split("\\s+")[5]);
		//	RT = s.spectrumName.split("\\s+")[5];
		//	RT = RT.substring(2, RT.length()-1);
		//	s.rt = Double.parseDouble(RT);	
		}
		if(s.spectrumName.contains("PROTEIN")){
			s.protein = s.spectrumName.split("PROTEIN: ")[1];
			//System.out.println("proteins: " + s.protein);
		}
		
		StringBuffer buff = new StringBuffer();
		buff.append("Name: " + s.peptide +"/" + s.charge+"\n");
		buff.append("MW: " + s.parentMass*s.charge +"\n");
		buff.append("Comment: ");
		buff.append("Fullname="+s.peptide+"/"+s.charge + " ");
		buff.append("Protein="+"\""+s.protein+"\""+" ");
		buff.append("RetentionTime="+"\""+RT+"\""+" ");
		buff.append("\n");
		buff.append("NumPeaks: " + s.getPeaks().size() +"\n");
		for(int i = 0; i < s.getPeak().size(); i++){
			Peak p = s.getPeak().get(i);
			buff.append(p.getMass() + "\t" + p.getIntensity()+"\t" + "\"? 1/1 0.0\""+"\n");
		}
		buff.append("\n");
		return buff.toString();
	}
	
	public static void testGenerateMSPLib(){
		String libFile = "../mixture_linked/Human_QTOF5600_tppSpecST_consensus_withRT.mgf";
		SpectrumWriter writer = new SpectrumWriter();
		SpectrumLib lib = new SpectrumLib(libFile, "MGF");
		writer.printSpectraToFile(lib.getAllSpectrums(), "..//mixture_linked//test_msp_out/test_HumanTOF5600.msp");
	}
	
	public static void main(String[] args){
		testGenerateMSPLib();
	}
	
	
	
}
