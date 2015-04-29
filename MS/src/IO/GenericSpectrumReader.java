package IO;

import org.Spectrums.Spectrum;


public class GenericSpectrumReader implements SpectrumReader{
	
	private String format;
	private String spectrumFile;
	private SpectrumReader reader;
	
	public GenericSpectrumReader(String spectrumFile, String format){
		this.format = format.toLowerCase();
		this.spectrumFile = spectrumFile;
		if(this.format.equals("mgf")){
			this.reader = new MGFReader(spectrumFile);
		}else if(this.format.equals("mzxml")){
			this.reader = new MZXMLReader(spectrumFile);
		}else if(this.format.equals("sptxt")){
			this.reader = null;
		}else if(this.format.equals("msp")){
			this.reader = null;
		}else if(this.format.equals("ms2")){
			this.reader = new MS2Reader(spectrumFile);
		}else{
			throw new IllegalArgumentException("Spectrum File format not supported: " + spectrumFile);
		}
	}
	
	public GenericSpectrumReader(String spectrumFile){
		this(spectrumFile, Utils.FileIOUtils.getFileExtension(spectrumFile));
	}
	
	
	@Override
	public boolean hasNext() {
		return reader.hasNext();
	}

	@Override
	public Spectrum next() {
		return reader.next();
	}

	@Override
	public void remove() {
		this.reader.remove();
	}

	@Override
	public Spectrum readSpectrumByIndex(int index) {
		return reader.readSpectrumByIndex(index);
	}
	
	public static void testSpectrumReader(){
		//String file = "..//mixture_linked//yeast_data//klc_010908p_yeast-digest.mzXML";
		String file = "..//mixture_linked//yeast_data//klc_010908p_yeast-digest.mgf";
		GenericSpectrumReader reader = new GenericSpectrumReader(file);
		while(reader.hasNext()){
			Spectrum s = reader.next();
			System.out.println("Spectrum " + s.scanNumber + " has peaks: " + s.getPeak().size());
		}
		int[] indices = new int[]{1, 10, 20, 30, 40, 50, 60, 70, 80};
		for(int i = 0; i < indices.length; i++){
			Spectrum s = reader.readSpectrumByIndex(indices[i]);
			System.out.println("Spectrum " + s.scanNumber + " has peaks: " + s.getPeak().size());
		}
		
		int[] indices2 = new int[]{1, 30, 20, 50, 70, 80, 40, 90, 180};
		for(int i = 0; i < indices.length; i++){
			Spectrum s = reader.readSpectrumByIndex(indices2[i]);
			System.out.println("Spectrum " + s.scanNumber + " has peaks: " + s.getPeak().size());
		}
	}
	
	public static void main(String[] args){
		testSpectrumReader();
	}
	
}
