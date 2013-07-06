package mixdb;

import java.util.Iterator;
import java.util.TreeMap;

import org.Spectrums.LargeHashMap;
import org.Spectrums.LargeSpectrumLibIterator;
import org.Spectrums.Spectrum;

/**
 * Read  a spectrum file in increasing parentmass manner
 * @author Jian Wang
 *
 */
public class SortedSpectrumReader implements Iterator{
	private String spectrumFile;
	private LargeHashMap specMap;
	private TreeMap<Double, Integer> indexMap;
	private Iterator<Integer> indexIter;
	public SortedSpectrumReader(String file){
		Iterator it = new LargeSpectrumLibIterator(file);
		this.spectrumFile = file;
		this.specMap = new LargeHashMap(spectrumFile+".spec");
		this.indexMap = new TreeMap<Double, Integer>();
		int currentIndex  = 1;
		while(it.hasNext()){
			Spectrum s = (Spectrum)it.next();
			if(s.scanNumber <= 0){
				s.scanNumber = currentIndex;
			}
			specMap.put(currentIndex, s);
			indexMap.put(s.parentMass+Math.random()*0.00001, currentIndex++);
		}
		specMap.finalize();
		this.indexIter = indexMap.values().iterator();
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
