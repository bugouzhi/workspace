package Utils;

import org.Spectrums.LabelledPeak;
import org.Spectrums.Peak;
import org.Spectrums.PeakIntensityComparator;
import org.Spectrums.PeakMassComparator;
import org.Spectrums.Peptide;
import org.Spectrums.SimpleMatchingGraph;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumLib;
import org.Spectrums.TheoreticalSpectrum;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * Contain utility methods to filter spectrum and methods
 * to compute other signal and nosie statistics for a spectrum
 * @author Jian
 *
 */
public class SpectrumFilter {
	public static int MODE_GLOBAL = 0;
	public static int MODE_LOCAL = 1;
	private int filteringMode = 1;
	private double lowFract = 0.5;
	private int topN = 8;
	private double windowWidth = 25.0;
	private double minSignalR = 2.5;
	private boolean globalFilter = true;
	private boolean localFilter = false;
	
	private Spectrum s;
	private List<Peak> signals;
	private List<Peak> noises;
	
	public boolean DEBUG=false;
	
	public int getFilteringMode() {
		return filteringMode;
	}

	public void setFilteringMode(int filteringMode) {
		this.filteringMode = filteringMode;
	}

	public double getLowFract() {
		return lowFract;
	}

	public void setLowFract(double lowFract) {
		this.lowFract = lowFract;
	}

	public int getTopN() {
		return topN;
	}

	public void setTopN(int topN) {
		this.topN = topN;
	}

	public double getWindowWidth() {
		return windowWidth;
	}

	public void setWindowWidth(double windowWidth) {
		this.windowWidth = windowWidth;
	}

	public double getMinSignalR() {
		return minSignalR;
	}

	public void setMinSignalR(double minSignalR) {
		this.minSignalR = minSignalR;
	}

	public List<Peak> getSignals() {
		return signals;
	}

	public void setSignals(List<Peak> signals) {
		this.signals = signals;
	}
	
	

	public boolean useGlobalFilter() {
		return globalFilter;
	}

	public void setGlobalFilter(boolean globalFilter) {
		this.globalFilter = globalFilter;
	}

	public boolean useLocalFilter() {
		return localFilter;
	}

	public void setLocalFilter(boolean localFilter) {
		this.localFilter = localFilter;
	}

	/**
	 * Define noise peaks as peaks with lowFract of lowest intesnity
	 * @param lowFract
	 * @return set of peaks designated as noises peaks
	 */
	public Set<Peak> splitByIntensity(Spectrum s){
		int numPeaks = s.getPeak().size();
		//System.out.println("Number of peaks " + numPeaks);
		int ind = (int)(lowFract*numPeaks);
		List<Peak> allPeaks = new ArrayList<Peak>();
		allPeaks.addAll(s.getPeak());
		Collections.sort(allPeaks, PeakIntensityComparator.comparator);
		List<Peak> signals = new ArrayList<Peak>();
		Set<Peak> noises = new HashSet<Peak>();
		for(int i = 0; i <= ind; i++){
			noises.add(allPeaks.get(i));
		}
		for(int i = ind+1; i < numPeaks; i++){
			signals.add(allPeaks.get(i));
		}
		//System.out.println("middle point is: " + ind);
		//System.out.println("signal size: " + signals.size() + "\tnoise size:\t" + noises.size());
		return noises;
	}
	
	/**
	 * Define noise as peak ranks in a local windows
	 * @param s
	 * @return set of peaks designated as noise peaks
	 */
	public Set<Peak> splitByLocalRank(Spectrum s){
		int current = 0, left = 0, right = 0;
		Set<Peak> noises = new HashSet();
		TreeSet<Peak> neighs = new TreeSet(PeakIntensityComparator.comparator);
		List<Peak> pList = new ArrayList<Peak>();
		pList.addAll(s.getPeaks());
		while(current < pList.size()){
			Peak p = pList.get(current);
			for(int i = left; i < current; i++){
				Peak smaller = pList.get(i);
				if(p.getMass() - smaller.getMass() > this.windowWidth){
					neighs.remove(smaller);
					left = i;
				}
			}
			for(int i = right+1;  i < pList.size(); i++){
				Peak bigger = pList.get(i);
				bigger.setIntensity(bigger.getIntensity()*(1+i*0.000001));   //creat some ordering for same intensity peaks
				if(bigger.getMass() - p.getMass()  <= this.windowWidth){
					neighs.add(bigger);
				}else{
					right = i-1;
					break;
				}
			}
			//System.out.println("current list size: " + neighs.size());
			Iterator<Peak> it = neighs.descendingIterator();
			Peak prev = null;
			for(int i = 0; it.hasNext();i++){
//				//System.out.println("i is: " +i);
				Peak p2 = it.next();
				//System.out.println(p2);
//				//System.out.println("peak " + p2);
				if(p == p2){
					break;
				}
				//if(prev != null && (prev.getIntensity() - p2.getIntensity() > 0.01)){
				//	i++;
				//}
				prev = p2;
				if(i > this.topN){
					noises.add(p);
					break;
//					System.out.println("remove " + p);
				}
			}
			//System.out.println();
			current++;
		}
	
		//System.out.println("signal size: " + pList.size() + "\tnoise size:\t" + noises.size());
		return noises;
	}
	
	public void filterSpectrum(Spectrum s, int mode){
		this.splitSpectrum(s);
		s.getPeak().clear();
		s.getPeak().addAll(this.signals);
		Collections.sort(s.getPeak(), PeakMassComparator.comparator);
		s.setPeaks(this.signals);
	}
	
	private void splitSpectrum(Spectrum s){
		Set<Peak> noises=new HashSet<Peak>();
		if(this.globalFilter){
			Set<Peak> noises1 = splitByIntensity(s);
			noises.addAll(noises1);
		}
		if(this.localFilter){
			Set<Peak> noises2 = splitByLocalRank(s);
			noises.addAll(noises2);
		}
		List<Peak> peaks  = new ArrayList<Peak>();
		peaks.addAll(s.getPeak());
		peaks.removeAll(noises);
		this.noises.addAll(noises);
		//System.out.println("size of noise " + this.noises.size());
		Collections.sort(this.noises, PeakMassComparator.comparator);
		this.signals = peaks;
		Collections.sort(this.signals, PeakMassComparator.comparator);
	}
	
	public double[] medianStat(Spectrum s){
		this.signals = new ArrayList<Peak>();
		this.noises = new ArrayList<Peak>();
		this.splitSpectrum(s);
		int index1 = (int)(this.signals.size()*0.5);
		int index2 = (int)(this.noises.size()*0.5);
		double medianNoise = 0.01, medianSignal=0.0;
		if(this.signals.size() > 0){
			medianSignal =this.signals.get(index1).getIntensity();
		}
		if(this.noises.size() > 0){
			medianNoise =this.noises.get(index2).getIntensity();
		}
		//System.out.println(this.signals.get(0).getIntensity() + "\t" + this.signals.get(this.signals.size()-1).getIntensity());
		//System.out.println(this.noises.get(0).getIntensity() + "\t" + this.noises.get(this.noises.size()-1).getIntensity());
		return new double[]{medianSignal, medianNoise};		
	}
	
	public  double[] computeSNR(Spectrum s){
		List<Peak>[] pList;
		int numSignals = 0;
		double[] SNRs = new double[s.getPeak().size()];
		double[] stat = medianStat(s);
		double medianSig = stat[0];
		double medianNoise = stat[1];
		
		//System.out.println("median noise: " + medianNoise);
		for(int i = 0; i < s.getPeak().size(); i++){
			double R = s.getPeak().get(i).getIntensity() / medianNoise;
			SNRs[i] = R;
			if(R > this.minSignalR){
				numSignals++;
			}
		}
		return new double[]{numSignals, medianNoise, medianSig/medianNoise};
	}
	
	
	public double[] computeAnnotatedSNR(Spectrum s, double tolerance){
		Peptide pep = new Peptide(s.peptide, 2);//s.charge);
		String[] prefixIons = {"b", "b-H20", "b-NH3", "a"};//, "b(iso)"};//"a", "a-H20", "a-NH3"};
		String[] suffixIons = {"y", "y-H20", "y-NH3"};//, "y(iso)"};
		TheoreticalSpectrum.prefixIons = prefixIons;
		TheoreticalSpectrum.suffixIons = suffixIons;
		TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
		double[] stats = new double[]{};
		double[] stat = medianStat(s);
		TheoreticalSpectrum.addIsotopicPeaks(t, 1);
		s = new Spectrum(s); //make a copy of the spectrum because there is some processing below that will alter the spectrum
		s.removePrecursors(0.1);
		SimpleMatchingGraph g = t.getMatchGraph(s, tolerance);
		List<Peak> annotated = new ArrayList();
		List<Peak> unAnnotated = new ArrayList();
		List<Peak> allPeaks = new ArrayList();
		allPeaks.addAll(s.getPeak());
		if(allPeaks.size() == 0){
			return stats;
		}
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				unAnnotated.add(p);
			}else{
				annotated.add(p);
			}
		}
		double medianNoise = medianStat(s)[1];
		int annotatedSig = 0, annotatedNoise = 0, unAnnotatedSig = 0, unAnnotatedNoise = 0;
		for(Iterator<Peak> it = annotated.iterator(); it.hasNext();){
			if(it.next().getIntensity()/medianNoise > this.minSignalR){
				annotatedSig++;
			}else{
				annotatedNoise++;
			}
		}
		for(Iterator<Peak> it = unAnnotated.iterator(); it.hasNext();){
			if(it.next().getIntensity()/medianNoise > this.minSignalR){
				unAnnotatedSig++;
			}else{
				unAnnotatedNoise++;
			}
		}
		double fractSigAnnotated = annotatedSig/(annotatedSig+unAnnotatedSig+0.01);
		double fractAnnotatedSig = annotatedSig / (annotatedSig + annotatedNoise+0.01);
		return new double[]{medianNoise, stat[0]/medianNoise, annotatedSig, fractSigAnnotated, fractAnnotatedSig, annotatedNoise, unAnnotatedNoise};
	}
	
	public static void testSpectrumFilter(){
		String spectrumFile = "../mixture_linked/test.mgf";
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		List<Spectrum> specList = lib.getAllSpectrums();
		SpectrumFilter filter = new SpectrumFilter();
		double[] numSignals = new double[specList.size()];
		double[] SNRs = new double[numSignals.length];
		for(int i = 0; i < specList.size(); i++){
			filter.setGlobalFilter(true);
			//filter.setLocalFilter(true);
			Spectrum s = specList.get(i);
			s.mergePeaks(s, 0.05);
			double[] stat = filter.computeAnnotatedSNR(s, 0.5);
			//double[] stat = filter.computeAnnotatedSNR(s, SpectrumFilter.MODE_GLOBAL, 0.05);
			System.out.println(s.spectrumName + "\t" + s.peptide + "\t" + s.parentMass + "\t" + s.charge + "\t" + stat[0] + "\t" + stat[1] + "\t" + stat[2]);
			numSignals[i] = stat[2];
			SNRs[i] = stat[1];
		}
		double[] signalBins = new double[]{0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,Integer.MAX_VALUE};
		double[] SNRBins = new double[]{0.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,10.0,Double.MAX_VALUE};
		int[] distr1 = ArrayUtils.hist1D(numSignals, signalBins);
		int[] distr2 = ArrayUtils.hist1D(SNRs, SNRBins);
		System.out.println(ArrayUtils.displayHist1D(distr1, signalBins));
		System.out.println();
		System.out.println(ArrayUtils.displayHist1D(distr2, SNRBins));
	}
	
	public static void main(String args[]){
		testSpectrumFilter();
	}
	
}
