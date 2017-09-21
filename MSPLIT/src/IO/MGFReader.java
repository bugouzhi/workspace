package IO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

import Spectrum.Spectrum;
/**
 * Reader of mgf file
 * @author Jian
 *
 */
public class MGFReader implements SpectrumReader{

	private String spectrumFile;
	private Spectrum nextSpectrum;
	private BufferedReader buff;
	private boolean isProceed=true;
	private boolean hasNext = false;
	private int beginInd = 1;
	
	//where does index begin for mgf file
	public int getBeginInd() {
		return beginInd;
	}

	public void setBeginInd(int beginInd) {
		this.beginInd = beginInd;
		reStart();
	}

	private int count = beginInd;
	public MGFReader(String filename){
		this.spectrumFile = filename;
		try{
			BufferedReader buff = new BufferedReader(new FileReader(this.spectrumFile));
			this.buff = buff;
			this.nextSpectrum = new Spectrum();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public List<Spectrum> readAllSpectra(){
		List<Spectrum> specList = new ArrayList<Spectrum>();
		int counter = 0;
		while(this.hasNext()){
			Spectrum s = (Spectrum)this.next();
			s.scanNumber = counter++;
			specList.add(s);
		}
		System.out.println("read in total spectra: " + specList.size());
		return specList;
	}
	
	@Override
	public boolean hasNext() {
		if(!isProceed){
			return hasNext;
		}
		if(nextSpectrum.readSpectrumFromMGF(buff)){
		//if(nextSpectrum.readSpectrumFromMS2(buff)){
			this.hasNext = true;
			isProceed = false;
			
			return true;
		}else{
			this.hasNext = false;
			this.nextSpectrum = new Spectrum();
			isProceed = false;
			return false;
		}
	}

	@Override
	public Spectrum next() {
		if(this.hasNext()){
			isProceed = true;
			Spectrum ret = this.nextSpectrum;
			this.nextSpectrum = new Spectrum();
			ret.scanNumber = this.count++;
			return ret;             //not very good way to do this
		}else{
			throw new NoSuchElementException();
		}
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}
	
	//restart the iterator
	public void reStart(){
		try{
			this.buff.close();
			BufferedReader buff = new BufferedReader(new FileReader(this.spectrumFile));
			this.buff = buff;
			this.nextSpectrum = new Spectrum();
			this.isProceed = true;
			this.hasNext = false;
			this.count = this.beginInd;
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}

	@Override
	public Spectrum readSpectrumByIndex(int index) {
		if(count == index){
			return next();
		}else if(count > index){
			reStart();
			return readSpectrumByIndex(index);
		}else if(count < index){
			while(count < index && hasNext()){
				next();
			}
			if(!hasNext()){
				throw new IllegalArgumentException("Spectrum index/scan not valide");
			}
			return this.nextSpectrum;
		}
		return null;
	}

	private Exception IllegalArgumentException(String string) {
		// TODO Auto-generated method stub
		return null;
	}

}
