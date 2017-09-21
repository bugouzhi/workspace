package IO;

import java.util.Iterator;
import java.util.TreeMap;

import Spectrum.Spectrum;
import Utils.LargeHashMap;


/**
 * A genreal interface to read  a spectrum file in increasing parentmass manner. 
 * The specific parser is called depending on the file type as judge by the extension
 * This is used during searches, so candidate spectra matches to query only need
 * to be loaded from storage once
 * @author Jian Wang
 *
 */
public class SortedSpectrumReader implements Iterator{
	private String spectrumFile;
	private LargeHashMap specMap;
	private TreeMap<Double, Integer> indexMap;
	private Iterator<Integer> indexIter;
	private double minParentMass = 0.0;
	
	public SortedSpectrumReader(String file){
		Iterator it = new JmzReaderAdaptor(file);
		this.spectrumFile = file;
		this.specMap = new LargeHashMap(spectrumFile+".spec");
		this.indexMap = new TreeMap<Double, Integer>();
		int currentIndex  = 1;
		System.out.println("Generateding sorted iterator");
		while(it.hasNext()){
			Spectrum s = (Spectrum)it.next();
			if(s.scanNumber <= 0){
				s.scanNumber = currentIndex;
			}
			specMap.put(currentIndex, s);
			//System.out.println("Putting " + s.spectrumName + "\t" + s.parentMass); 
			if(s.parentMass > this.minParentMass) {  //skipping ms1 spectra
				indexMap.put(s.parentMass+Math.random()*0.00001, currentIndex++);
			}
			if(currentIndex % 10000 == 0){
				System.out.println("Processed spectrum " + currentIndex);
			}
		}
		specMap.finalizeObjectFile();
		this.indexIter = indexMap.values().iterator();
		System.out.println("Done with sorted iterator");
	}
	@Override
	public boolean hasNext() {
		return indexIter.hasNext();
	}
	@Override
	public Object next() {
		int index = this.indexIter.next();
		return (Spectrum)specMap.get(index);
	}
	@Override
	public void remove() {
		System.err.println("warning method not implemented yet, it does no modification now");
		// TODO Auto-generated method stub
		
	}
	
	
	public double getMinParentMass() {
		return minParentMass;
	}
	public void setMinParentMass(double minParentMass) {
		this.minParentMass = minParentMass;
	}
	
	public static void testSortedReader(){
		String spectrumFile = "../mixture_linked/yeast_simulatedmixture_samecharge_alpha1.0.mgf";
		SortedSpectrumReader it = new SortedSpectrumReader(spectrumFile);
		while(it.hasNext()){
			Spectrum s = (Spectrum)it.next();
			System.out.println("spectrum: " + s.scanNumber + "\t" + s.parentMass + "\t" + s.charge + "\t" + s.getPeak().size());
		}
	}
	
	public static void main(String[] args){
		testSortedReader();
	}

}
