package MSPLIT.SpectrumIO;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import MSPLIT.Spectrum;

/**
 * General class for parsing a spectra file
 * 
 * @author bugouzhi
 *
 */
public abstract class SpectrumParser implements Iterator<Spectrum>{	
	public String getSpectrumFile() {
		return spectrumFile;
	}

	public void setSpectrumFile(String spectrumFile) {
		this.spectrumFile = spectrumFile;
	}

	protected BufferedReader getBuff() {
		return buff;
	}

	public void setBuff(BufferedReader buff) {
		this.buff = buff;
	}

		private String spectrumFile;
		private BufferedReader buff;
		private Spectrum nextSpectrum;
		
		public SpectrumParser(String filename){
			this.spectrumFile = filename;
			try{
				BufferedReader buff = new BufferedReader(new FileReader(this.spectrumFile));
				this.buff = buff;
				this.nextSpectrum = this.getNextSpectrum();
			}catch(IOException ioe){
				System.out.println(ioe.getMessage());
				ioe.printStackTrace();
			}
		}
		
		public List<Spectrum> readAllSpectra(){
			List<Spectrum> list = new ArrayList<Spectrum>();
			while(this.hasNext()){
				Spectrum s = this.next();
				list.add(s);
			}
			return list;
		}
		@Override
		public boolean hasNext() {
			return nextSpectrum != null;
		}
		
		@Override
		public Spectrum next() {
			Spectrum ret = nextSpectrum;
			this.nextSpectrum = getNextSpectrum(); 
			return ret;
		}
		
		protected abstract Spectrum getNextSpectrum();
		
		@Override
		public void remove() {
			//not supporting  remove in this case
			// TODO Auto-generated method stub
			
		}
		
}
