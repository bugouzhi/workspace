package IO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import Spectrum.Spectrum;

/**
 * When spectrum lib is large we avoid loading all of them into memory at once
 * instead we read from file one spectrum at a time, the user process through the
 * spectrum through a iterator interface
 * @author Jian Wang
 *
 */
public class LargeSpectrumLibIterator<T> implements Iterator<T>{
	private String spectrumFile;
	private Spectrum nextSpectrum;
	private BufferedReader buff;
	private boolean isProceed=true;
	private boolean hasNext = false;
	private int count = 1;
	public LargeSpectrumLibIterator(String filename){
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
	public T next() {
		if(this.hasNext()){
			isProceed = true;
			Spectrum ret = this.nextSpectrum;
			ret.specIndex = count++;
			this.nextSpectrum = new Spectrum();
			return (T)ret;             //not very good way to do this
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
			BufferedReader buff = new BufferedReader(new FileReader(this.spectrumFile));
			this.buff = buff;
			this.nextSpectrum = new Spectrum();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}

}
