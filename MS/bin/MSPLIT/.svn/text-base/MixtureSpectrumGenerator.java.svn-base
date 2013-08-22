package MSPLIT;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Vector;


/**
 * Simulate generation of mixture spectrum
 * @author bugouzhi
 *
 */
public class MixtureSpectrumGenerator{
	private double minScale=1;
	private double maxScale=1;
	private double minSim = 0.0;   
	private double maxSim = 1.0;
	private double precursorMassDiff=3.0;
	private SpectrumComparator comp = CosineSpectrumComparator.CosineComparator;
	
	private List<AnnotatedSpectrum> specList; //internal list of spectra to generate mixture spectra
	
	public MixtureSpectrumGenerator(List<AnnotatedSpectrum> specList){
		this.specList = specList;
	}
	
	public MixtureSpectrumGenerator(Iterator<AnnotatedSpectrum> specIter){
		this.specList = new ArrayList<AnnotatedSpectrum>();
		while(specIter.hasNext()){
			this.specList.add(specIter.next());
		}
	}
	
	public MixtureSpectrumGenerator(Iterator<AnnotatedSpectrum> specIter,
			double scale1, double scale2, double minSim, double maxSim, double precursorDiff){
		this(specIter);
		this.minScale = scale1;
		this.maxScale = scale2;
		this.minSim = minSim;
		this.maxSim = maxSim;
		this.precursorMassDiff = precursorDiff;
	}
	
	/**
	 * randomly select spectrum from the list of spectra and
	 * generate mixture spectrum
	 * @param count
	 * @return
	 */
	public List<Spectrum> generateRandomMixtures(int count){
		return generateRandomMixtures(count, minScale, maxScale);
	}
	
	public List<Spectrum> generateRandomMixtures(int count, double minScale, double maxScale){
		int index1 = 0, index2 = 0, currentCount = 0;
		int size = this.specList.size()-1;
		double alpha = 1.0;
		this.minScale = minScale;
		this.maxScale = maxScale;
		List<Spectrum> mixtures = new ArrayList<Spectrum>();
		while(currentCount < count){
			index1 = (int)(Math.random()*size);
			index2 = (int)(Math.random()*size);
			if(index1 != index2){
				AnnotatedSpectrum s1 = this.specList.get(index1);
				AnnotatedSpectrum s2 = this.specList.get(index2);
				//System.out.println("checking");
				if(checkValidMixturePair(s1, s2)){
					alpha = this.generateAlpha();
					mixtures.add(new AnnotatedMixtureSpectrum(s1, s2, 1, alpha));
					currentCount++;
				}
			}
		}	
		return mixtures;
	}
	
	/**
	 * Iterate through the spectra list and generate mixtures until the 
	 * specified number of mixture spectra are generated
	 * @param count
	 * @return
	 */
	public List<Spectrum> generateMixtures(int count){
		return generateMixtures(count, minScale, maxScale);
	}
	
	public List<Spectrum> generateMixtures(int count, double minScale, double maxScale){
		int index1 = 0, index2 = 0, currentCount = 0;
		int size = this.specList.size()-1;
		double alpha = 1.0;
		this.minScale = minScale;
		this.maxScale = maxScale;
		List<Spectrum> mixtures = new ArrayList<Spectrum>();
		while(currentCount < count){
			if(index2 == size-1){
				index1++;
				index2 = index1+1;
			}else{
				index2++;
			}
			
			if(index1 != index2){
				AnnotatedSpectrum s1 = this.specList.get(index1);
				AnnotatedSpectrum s2 = this.specList.get(index2);
				if(checkValidMixturePair(s1, s2)){
					alpha = generateAlpha();
					mixtures.add(new AnnotatedMixtureSpectrum(s1, s2, 1.0, alpha));
					System.out.println(currentCount);
					currentCount++;
				}
			}
		}	
		return mixtures;
	}
	
	private boolean checkValidMixturePair(Spectrum s1, Spectrum s2){
		//System.out.println("checking valid pair");
		double sim = this.comp.compare(s1, s2);
		return sim > this.minSim && sim < this.maxSim 
			&& Math.abs(s1.getParentMass() - s2.getParentMass()) < this.precursorMassDiff;
	}
	
	private double generateAlpha(){
		return minScale + (maxScale-minScale)*Math.random();
	}
}
