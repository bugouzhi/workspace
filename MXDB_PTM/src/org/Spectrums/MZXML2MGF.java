package org.Spectrums;

import java.io.File;
import java.util.List;

/**
 * Convert mzxml file to mgf
 * @author Jian Wang
 *
 */
public class MZXML2MGF {
	public static void convertToMGF(String spectrumFile){
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		String mgfname = spectrumFile.substring(0, spectrumFile.lastIndexOf('.'));
		mgfname = mgfname+".mgf";
		System.out.println("In spectrumFile: " + spectrumFile +"\tout spectrumfile: " + mgfname);
		SpectrumUtil.printSpectraToFile(mgfname, specList);
	}
	
	public static void convertToMGFS(String spectrumDir){
		File f = new File(spectrumDir);
		String[] files=null;
		System.out.println("directory is: " + spectrumDir + "\t" + f.isDirectory());
		if(f.isDirectory()){
			files=f.list();
			for(int i = 0; i < files.length; i++){
				if(files[i].contains("mzXML")){
					convertToMGF(spectrumDir+"/"+files[i]);
				}
			}

		}
	}
	
	public static void main(String[] args){
		convertToMGF("..\\mixture_linked\\msdata\\linked_peptide_library\\ACG_disulfide_library\\08_30_2012\\PepLib1_TCEP_300ng.mzXML");
		//convertToMGFS("../mixture_linked/msdata/Ubiquitynlation/Yeast/mock/");
	}
}
