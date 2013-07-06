package MSPLIT;
import java.util.Iterator;
import java.util.List;
import java.io.Serializable;
import java.util.ArrayList;

/**
 * An implementation of spectrum back by a List
 * @author jian wang
 *
 */
public class SimpleSpectrum implements Spectrum, Serializable{
	public static final long serialVersionUID = 1L;
	//the peaks in this spectrum, we should always keep this list sorted
	//according to the mass
	protected List <Peak> peaks;
	protected Peak[] peaksArry;
	private String spectrumName = "";
	private double parentMass;
	private int charge;
	
	
	public SimpleSpectrum(){
		peaks = new ArrayList<Peak>();
		spectrumName = "DefaultSpectrum";
		parentMass = 100.0;
		charge = 1;
	}
	
	public SimpleSpectrum(String spectName, double parentMass, int charge, List<Peak> peakList){
		this.setSpectrumName(spectName);
		this.setParentMass(parentMass);
		this.setCharge(charge);
		//this.peaks = peakList;
//		this.peaksArry = new Peak[peakList.size()];
//		for(int i = 0; i < peakList.size(); i++){
//			peaksArry[i] = peakList.get(i);
//		}
	}
	
	//copy constructor
	public SimpleSpectrum(Spectrum s){
		this(s.getSpectrumName(), s.getParentMass(), s.getCharge(), new ArrayList<Peak>());
		for(int i = 0; i < s.numOfPeaks(); i++){
			this.addPeak(s.getPeak(i));               //shared peaks may not be desirable
		}
	}

	@Override
	public void addPeak(Peak p) {
		if(p == null){
			throw new IllegalArgumentException("peak cannot be null");
		}
		if(this.peaks == null){
			this.peaks = new ArrayList<Peak>();
		}
		for(int i = 0, size = peaks.size(); i <= size; i++){
			if(i == size){
				peaks.add(p);
				break;
			}
			Peak current = peaks.get(i);
			if(current.getMass() >=  p.getMass()){
				peaks.add(i, p);
			}
		}
		
	}
	
	/**
	 * return a list view of all the peaks
	 * @return
	 */
	public List<Peak> getAllPeaks() {
		//List<Peak> copy = new ArrayList<Peak>();
		//scopy.addAll(peaks);
		return peaks;
	}
	
	@Override
	public Peak getPeak(int position) {
		return peaks.get(position);
		//return this.peaksArry[position];
	}
	
	@Override
	public int numOfPeaks() {
		return peaks.size();
	}
	@Override
	public Iterator<Peak> peakIterator() {
		return peaks.iterator();
	}
	
	@Override
	public Peak removePeak(int i) {
		return peaks.remove(i);
	}
	@Override
	public boolean removePeak(Peak p) {
		return peaks.remove(p);
	}

	public String getSpectrumName() {
		return spectrumName;
	}

	public void setSpectrumName(String spectrumName) {
		this.spectrumName = spectrumName;
	}

	public double getParentMass() {
		return parentMass;
	}

	public void setParentMass(double parentMass) {
		this.parentMass = parentMass;
	}

	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		if(charge < 0){
			throw new IllegalArgumentException("charge cannot be less than zero");
		}
		this.charge = charge;
	}
	
	public void setPeaks(List<Peak> peakList) {
		this.peaks = peakList;
	}
	
}
