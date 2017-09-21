package IO;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import Spectrum.Peak;
import Spectrum.PeakMassComparator;

import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;
//import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzml_wrapper.MzMlWrapper;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
import Spectrum.Spectrum;

/**
 * This class provide a adapting interface to use the jmzReader interface, converting spectrum
 * object from the jmzReader framework to our own spectrum data object. Since the jmzReader supports
 * several standard format (e.g. mzml, mzxml), we use their parser instead of the built-in one in msplit to extends the support
 * of query spectra format.  However for reading annotated spectrum from library (e.g. SEQ field etc..), the default interface of jzmReader
 * does not support these field yet since there are not standard terms/fields accross all formats.  The data are available, just the individual
 * application needs get these data depending on file format.  For now to load annotated files use the SpectrumReader interface.
 * @author Jian
 *
 */
public class JmzReaderAdaptor implements Iterator{
	private File file;
	private String format;
	private JMzReader reader;
	private Iterator it;
	private int count=1;
	private int minLevel = 2;
	
	/**
	 * 
	 * @param file
	 */
	public JmzReaderAdaptor(String file){
		this.file  = new File(file);
		this.format = Utils.FileIOUtils.getFileExtension(file);	
		this.format = this.format.toLowerCase();
		try{
			System.out.println("Spectrum file format: " + format);
			if(this.format.equals("mgf")){
				this.reader = new MgfFile(this.file);
			}else if(this.format.equals("mzxml")){
				this.reader = new MzXMLFile(this.file);
			}else if(this.format.equals("mzml")){
				this.reader = new MzMlWrapper(this.file);
			}else{
				throw new IllegalArgumentException("Unsupported spectrum format: " + this.format);
			}
			this.it = this.reader.getSpectrumIterator();
		}catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
	}

	@Override
	public boolean hasNext() {
		return this.it.hasNext();
	}

	@Override
	public Object next() {
			//Spectrum s = (Spectrum)this.it.next();
		if(this.it.hasNext()){
			uk.ac.ebi.pride.tools.jmzreader.model.Spectrum s = null;
			s = (uk.ac.ebi.pride.tools.jmzreader.model.Spectrum)this.it.next();
			return convert(s);
		}else{
			return null;
		}
	}	
	@Override
	public void remove() {
		this.it.remove();
		
	}
	
	public int getSpectrumCount(){
		return this.reader.getSpectraCount();
	}
	
	public Spectrum getSpectrumByIndex(int ind){
		try{
			return convert(this.reader.getSpectrumByIndex(ind));
		}catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
		return null;
	}
	
	public Spectrum getSpectrumById(String id){
		try{
			return convert(this.reader.getSpectrumById(id));
		}catch(JMzReaderException e){
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
		return null;
	}
	
	public Spectrum getSpectrumByScan(int scan){
		try{
			return convert(this.reader.getSpectrumById(""+scan));
		}catch(JMzReaderException e){
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
		return null;
	}
	
	private Spectrum convert(uk.ac.ebi.pride.tools.jmzreader.model.Spectrum s){
		//System.out.println(s.getId() + "\t" + s.getPrecursorMZ() + "\t" + s.getPrecursorCharge() +"\t" + s.getMsLevel());
		Spectrum converted = new Spectrum();
		converted.spectrumName = s.getId();
		//converted.specIndex = s.get;
		if(this.format.equals("mzxml")){
			converted.scanNumber = Integer.parseInt(s.getId());
		}
		if(this.format.equals("mgf")){
			converted.scanNumber = Integer.parseInt(((Ms2Query)s).getScan());
		}
		converted.specIndex = this.count++;
		converted.parentMass = s.getPrecursorMZ() == null ? 0 : s.getPrecursorMZ();
		converted.charge = s.getPrecursorCharge() == null ? 0 : s.getPrecursorCharge();
		Map<Double, Double> peakList = s.getPeakList();
		if(peakList==null){
			return converted;
		}
		//System.out.println("peaklist size: " + s.getPeakList().size());
		List<Peak> peaks = new ArrayList<Peak>();
		for(Iterator<Double> it = peakList.keySet().iterator(); it.hasNext();){
			double mz = it.next();
			peaks.add(new Peak(mz, peakList.get(mz)));
		}
		Collections.sort(peaks, PeakMassComparator.comparator);
		converted.setPeaks(peaks);
		return converted;
	}
	
	public static void testSpectrumReaderAdaptor(){
		String MgfFile = "C:\\Users\\Jian\\workspace\\mixture_linked\\MSPLIT_v1.0/test_readMGF/specs_ms_0.mgf";
		JmzReaderAdaptor parser = new JmzReaderAdaptor(MgfFile);
		System.out.println("Mgf file has spectra: " + parser.getSpectrumCount());
		
		String MzxmlFile = "..//mixture_linked//18487_REP3_40fmol_UPS1_1ug_Ecoli_NewStock2_IDA_1_subset.mzXML";
		//parser = new JmzReaderAdaptor(MzxmlFile);
		//System.out.println("Mzxml file has spectra: " + parser.getSpectrumCount());
		
		//String MzmlFile = "C:\\Users\\Jian\\workspace\\mixture_linked\\yeast_data/klc_010908p_yeast-digest.mzML";
		//parser = new JmzReaderAdaptor(MzmlFile);
		//System.out.println("Mzml file has spectra: " + parser.getSpectrumCount());
		int count= 0;
		while(parser.hasNext()){
			Spectrum s = (Spectrum)parser.next();
			System.out.println(s.spectrumName + "\t" + "\tindex\t" + s.specIndex + "\t" + s.parentMass + "\t" + s.charge + "\t" + s.getPeak().size() +"\t" + s.scanNumber);
			if(count > 10){
				break;
			}
			count++;
		}
	}
	
	public static void main(String[] args){
		testSpectrumReaderAdaptor();
	}
}
