package IO;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import Spectrum.Peak;
import Spectrum.PeakMassComparator;

//import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

import Spectrum.Spectrum;

public class MZMLReader implements SpectrumReader{
	File xmlFile;
	MzMLUnmarshaller unmarshaller;
	MzMLObjectIterator spectrumIterator;
	public MZMLReader(String file){
		this.xmlFile = new File("path/to/your/mzml/file");
		this.unmarshaller = new MzMLUnmarshaller(xmlFile);
		this.spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);
//		while (spectrumIterator.hasNext()){
//		  //read next spectrum from XML file
//		  Spectrum spectrum = spectrumIterator.next();
//		  //use it
//		  System.out.println("Spectrum ID: " + spectrum.getId());
//		}
	}
	
	
	public Spectrum readSpectrumByIndex(int index){
		Spectrum s = new Spectrum();
		return s;
	}

	@Override
	public boolean hasNext() {
		this.spectrumIterator.hasNext();
		return false;
	}

	@Override
	public Spectrum next() {
		if(hasNext()){
			uk.ac.ebi.jmzml.model.mzml.Spectrum s = (uk.ac.ebi.jmzml.model.mzml.Spectrum)this.spectrumIterator.next();
			return convert(s);
		}
		return null;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}
	
	private Spectrum convert(uk.ac.ebi.jmzml.model.mzml.Spectrum s){
		Spectrum converted = new Spectrum();
		converted.spectrumName = s.getId();
		converted.specIndex = s.getIndex();
		//converted.parentMass = s.
		//converted.specIndex = s.get;
		//converted.parentMass = s.getPrecursorMZ() == null ? 0 : s.getPrecursorMZ();
		//converted.charge = s.getPrecursorCharge() == null ? 0 : s.getPrecursorCharge();
		//Map<Double, Double> peakList = s.getPeakList();
		//List<Peak> peaks = new ArrayList<Peak>();
		//for(Iterator<Double> it = peakList.keySet().iterator(); it.hasNext();){
		//	double mz = it.next();
		//	peaks.add(new Peak(mz, peakList.get(mz)));
		//}
		//Collections.sort(peaks, PeakMassComparator.comparator);
		return converted;
	}
}
