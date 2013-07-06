package org.Spectrums;

import java.io.File;
import java.util.List;

import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.MzML;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.Scan;
import uk.ac.ebi.jmzml.model.mzml.ScanList;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;

/**
 * MZML reader for mass spectrum	
 * @author Jian Wang
 *
 */
public class MZMLReader {
	
	public static void testRunMZMLReader(){
		String spectrumFile = "../mixture_linked/msdata/linked_peptide_library/ACG_disulfide_library/Orbi_Elite/PepLib1_300ng_trp_Elite_CID_35rep.mzML"; 
		File xmlFile = new File(spectrumFile);
		MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(xmlFile);
		//MzML completeMzML = unmarshaller.unmarshall();
		MzMLObjectIterator<Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);
		double ms1time = 0.0;
		while (spectrumIterator.hasNext()){
		  //read next spectrum from XML file
		  Spectrum spectrum = spectrumIterator.next();
		  Spectrum ms1;
		  if(spectrum.getPrecursorList() == null){
			  ms1 = spectrum;
			  ms1time = getIonInjection(ms1);
			  continue;
		  }
		  double precursorInt = getPrecursorInt(spectrum);
		  double totalCurrent = Double.parseDouble(getParam(spectrum.getCvParam(), "total ion current"));
		  double ms2time = getIonInjection(spectrum);
		  System.out.println(spectrum.getId() + "\t" + precursorInt + "\t" + ms1time + "\t" + totalCurrent + "\t" + ms2time);
		}	
	}
	

	private static double getIonInjection(Spectrum spectrum){
		ScanList scanList = spectrum.getScanList();
		Scan scan = scanList.getScan().get(0);
		return Double.parseDouble(getParam(scan.getCvParam(), "ion injection time"));
	}
	
	private static double getPrecursorInt(Spectrum s){
		Precursor p = s.getPrecursorList().getPrecursor().get(0);
		String pInt = getParam(p.getSelectedIonList().getSelectedIon().get(0).getCvParam(), "peak intensity");
		return Double.parseDouble(pInt);
	}
	
	private static String getParam(List<CVParam> paramList, String paramName){
		  for(int i = 0; i < paramList.size(); i++){
			  CVParam param = paramList.get(i);
			  if(param.getName().equals(paramName)){
				  return param.getValue();
			  }
		  }
		 return "";
	}

	
	public static void main(String[] args){
		testRunMZMLReader();
	}
}
