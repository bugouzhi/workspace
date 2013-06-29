package mixdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.Spectrums.LabelledPeak;
import org.Spectrums.Mass;
import org.Spectrums.Peak;
import org.Spectrums.PeakMassComparator;
import org.Spectrums.Peptide;
import org.Spectrums.PeptideLite;
import org.Spectrums.Spectrum;
import org.Spectrums.TheoreticalSpectrum;
import sequences.FastaSequence;

/**
 * Factory to create theoretcial spectrum from single peptide
 * @author Jian Wang
 *
 */
public class TheoreticalSpectrumFactory {
	/**
	 * A type map for the standard ion types in the old system, can use to create
	 * theoretical spectrum from a peptide.
	 */
	public static Map<String, IonType> standardTypeMap = TheoreticalSpectrumFactory.createStandardIonTypeMap();
	public static IonTypeMapper standardIonMap = IonTypeMapper.createIonTypeMap(standardTypeMap.values());
	public static Map<String, Collection<IonType>> peptideMap = createPeptideMap(standardTypeMap);
	public static TheoreticalSpectrum getTheoreticalSpectrum(Peptide p){
		return new TheoreticalSpectrum(p);
	}
	
	public static TheoreticalSpectrum getTheoreticalSpectrum(PeptideLite p, FastaSequence protein, int charge){
		Peptide pep = new Peptide(protein.getSubsequence(p.getBeginInd(), p.getEndInd()+1));
		return new TheoreticalSpectrum(pep);
	}
	
	public static ArrayTheoreticalSpectrum getArrayTheoSpectrum(PeptideLite p, FastaSequence protein, int charge){
		Peptide pep = new Peptide(protein.getSubsequence(p.getBeginInd(), p.getEndInd()+1));
		return getArrayTheoSpectrum(pep, standardTypeMap, standardIonMap);
	}
	
	public static int getPeptideLength(Peptide p){
		return getPeptideLength(p.getPeptide(), p.getCharge());
	}
	
	public static int getPeptideLength(String p, int charge){
		if(charge == 2){
			if(p.length() <= 12){
				return SimplePeptideType.SHORT;
			}else{
				return SimplePeptideType.LONG;
			}
		}else if(charge == 3){
			if(p.length() <= 12){
				return SimplePeptideType.SHORT;
			}else{
				return SimplePeptideType.LONG;
			}
		}
		return 0;
	}
	
	public static Map<String, IonType> createStandardIonTypeMap(){
		int[] peptideLength = {SimplePeptideType.SHORT,SimplePeptideType.LONG};
		int[] peptideCharge = {2,3,4};
		String[] pTypes = TheoreticalSpectrum.prefixIons;
		String[] sTypes = TheoreticalSpectrum.suffixIons;
		Map<String, IonType> typeMap = new LinkedHashMap<String, IonType>();
		for(int l = 0; l < peptideLength.length; l++){
			for(int p = 0; p < peptideCharge.length; p++){
				PeptideType pepType = new SimplePeptideType(peptideCharge[p], peptideLength[l]);
				for(int c = 1; c < 5; c++){
					for(int t = 0; t < pTypes.length; t++){
						IonType ionType = new SimpleIonType(pTypes[t], Mass.getIonMod(pTypes[t]), c, SimpleIonType.PREFIX, pepType);
						typeMap.put(pTypes[t]+"@"+c+"@"+peptideCharge[p]+"@"+peptideLength[l], ionType);
					}
					
					for(int t = 0; t < sTypes.length; t++){
						IonType ionType = new SimpleIonType(sTypes[t], Mass.getIonMod(sTypes[t]), c, SimpleIonType.SUFFIX, pepType);
						typeMap.put(sTypes[t]+"@"+c+"@"+peptideCharge[p]+"@"+peptideLength[l], ionType);
					}
				}
			}
		}
		return typeMap;
	}
	
	public static Map<String, Collection<IonType>> createPeptideMap(Map<String, IonType> typeMap){
		Map<String, Collection<IonType>> peptideMap = new HashMap<String, Collection<IonType>>();
		for(Iterator<IonType> it = typeMap.values().iterator(); it.hasNext();){
			SimpleIonType type = (SimpleIonType)it.next();
			SimplePeptideType pType = (SimplePeptideType)type.getPType();
			String key = pType.toString();
			//System.out.println("key is: " + key);
			Collection<IonType> types;
			if(peptideMap.containsKey(key)){
				types = peptideMap.get(key);
			}else{
				types = new ArrayList<IonType>();
			}
			if(type.getCharge() <= pType.getPepCharge()){
				types.add(type);
			}
			peptideMap.put(key, types);
		}
		return peptideMap;
	}
	
	public static ArrayTheoreticalSpectrum getArrayTheoSpectrum(Peptide pep, Map<String, IonType> typeMap, IonTypeMapper ionMap){
		TheoreticalSpectrum th = new TheoreticalSpectrum(pep);
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] massIntList = new double[2][th.getPeak().size()];
		for(int i = 0; i < th.getPeak().size(); i++){
			LabelledPeak lp = (LabelledPeak)th.getPeak().get(i);
			String type = lp.getType()+"@"+lp.getCharge()+"@"+lp.getPep().getCharge()+"@"+getPeptideLength(pep);
			if(typeMap.containsKey(type)){
				int index = ionMap.getIndex(typeMap.get(type));
				massIntList[ArraySpectrum.MASS][i] = lp.getMass();
				massIntList[ArraySpectrum.INTENSITY][i] = index;
			}else{
				massIntList[ArraySpectrum.MASS][i] = lp.getMass();
				massIntList[ArraySpectrum.INTENSITY][i] = -1;
				System.err.println("warnining: ion-type " + type + " is not known in type map, treated as unknown type: -1");
			}
		}
		arry.setMassIntensityList(massIntList);
		System.out.println("spectrum has size: " + massIntList[0].length);
		arry.setPeptide(pep.getPeptide()+"."+pep.getCharge());
		arry.setCharge(pep.getCharge());
		arry.parentMass = pep.getParentmass();
		return arry;
	}
	
	public static ArrayTheoreticalSpectrum getArrayTheoSpectrum(PeptideLite pep, FastaSequence protein, Map<String, IonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		
		return arry;
	}
	

	public static Spectrum getTheoSpectrum(PeptideLite pep, FastaSequence protein, int charge, Map<String, IonType> typeMap, IonTypeMapper ionMap){
		String peptide = protein.getSubsequence(pep.getBeginInd(), pep.getEndInd()+1);
		return getTheoSpectrum(peptide, charge, typeMap, ionMap);
	}
	
	public static Spectrum getTheoSpectrum(String peptide, int charge, Map<String, IonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] baseMasses = computeBaseMass(peptide, new int[]{}, new double[]{});
		String[] leng = new String[]{"short", "long"};
		Collection<IonType> types = peptideMap.get(""+charge+"@"+leng[getPeptideLength(peptide, charge)]);
		//System.out.println("type length: " + types.size());
		List<Peak> peaks = new ArrayList<Peak>();
		TreeSet peaks2 = new TreeSet(PeakMassComparator.comparator);
		for(Iterator<IonType> it = types.iterator(); it.hasNext();){
			//for(int t = 0; t < 144; t++){
			SimpleIonType type = (SimpleIonType)it.next();
			//System.out.println("type is: " + type);
			int index = standardIonMap.getIndex(type);
			if(type.getDirection() == SimpleIonType.PREFIX){
				for(int i = 0; i < baseMasses[0].length; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*peakCharge)/peakCharge;
					Peak p = new Peak(mass, index);
					//peaks.add(p);
					//peaks2.add(p);
				}
			}		
			if(type.getDirection() == SimpleIonType.SUFFIX){
				for(int i = 0; i < baseMasses[1].length; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[1][i]+type.getOffset()+Mass.PROTON_MASS*peakCharge)/peakCharge;
					Peak p = new Peak(mass, index);
					//peaks.add(p);
					//peaks2.add(p);
				}
			}
		}
		//Collections.sort(peaks, PeakMassComparator.comparator);
		//arry.setPeaks(peaks);
		return arry;
	}
	
	
	public static List<Spectrum> getArrayTheoreticalSpectrum(List<PeptideLite> peps, FastaSequence seq, Spectrum s, int charge){
		List candidates = new ArrayList();
		Set<String> uniquePeps = new HashSet<String>(peps.size());
		for(int i = 0; i < peps.size(); i++){
			PeptideLite pep = peps.get(i);
			String pseq = seq.getSubsequence(pep.getBeginInd(), pep.getEndInd()+1);
			if(uniquePeps.contains(pseq)){  //check for redundancy
				continue;
			}
			uniquePeps.add(pseq);
			Peptide p = new Peptide(pseq+"."+charge);
			ArrayTheoreticalSpectrum th = (ArrayTheoreticalSpectrum)TheoreticalSpectrumFactory.getTheoSpectrumX(pep, seq, p.getCharge(), 
					TheoreticalSpectrumFactory.standardTypeMap, 
					TheoreticalSpectrumFactory.standardIonMap);
			
			candidates.add(th);
		}
		return candidates;
	}
	
	public static Spectrum getTheoSpectrumX(PeptideLite peplite, FastaSequence seq, int charge, Map<String, IonType> typeMap, IonTypeMapper ionMap){
		String pep=seq.getSubsequence(peplite.getBeginInd(), peplite.getEndInd()+1);
		ArrayTheoreticalSpectrum theo = (ArrayTheoreticalSpectrum)getTheoSpectrumX(pep, charge, typeMap, ionMap);
		peplite.setFastaseq(seq);
		theo.setPeplite(peplite);
		return theo;
	}
	
	//a fast array-base implementation
	public static Spectrum getTheoSpectrumX(String peptide, int charge, Map<String, IonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] baseMasses = computeBaseMass(peptide, new int[]{}, new double[]{}); 
		String[] leng = new String[]{"short", "long"};
		Collection<IonType> types = peptideMap.get(""+charge+"@"+leng[getPeptideLength(peptide, charge)]);
		double[][] massIntList = new double[2][baseMasses[0].length*types.size()];
		int currentIndex = 0;
		for(Iterator<IonType> it = types.iterator(); it.hasNext();){
			SimpleIonType type = (SimpleIonType)it.next();
			int index = standardIonMap.getIndex(type);
			if(type.getDirection() == SimpleIonType.PREFIX){
				for(int i = 0; i < baseMasses[0].length; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntList[0][currentIndex] = mass;
					massIntList[1][currentIndex++] = index;
				}
			}		
			if(type.getDirection() == SimpleIonType.SUFFIX){
				for(int i = 0; i < baseMasses[1].length; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[1][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntList[0][currentIndex] = mass;
					massIntList[1][currentIndex++] = index;
				}
			}
		}
		//Collections.sort(peaks, PeakMassComparator.comparator);
		//arry.setPeaks(peaks);
		//trying to sort by masses
		//double[][] sorted = sortByMass(massIntList);
		quicksort(massIntList[0], massIntList[1]);
//		for(int i = 0; i < massIntList[0].length; i++){
//			System.out.print(massIntList[0][i]+"\t");
//		}
//		System.out.println();
		Peptide p = new Peptide(peptide+"."+charge);
		arry.setMassIntensityList(massIntList);
		arry.peptide = peptide+"."+charge;
		arry.parentMass = p.getParentmass();
		arry.charge = charge;
		return arry;
	}
	
	public static void quicksort(double[] main, double[] index) {
	    quicksort(main, index, 0, index.length - 1);
	}

	// quicksort a[left] to a[right]
	public static void quicksort(double[] a, double[] index, int left, int right) {
	    if (right <= left) return;
	    int i = partition(a, index, left, right);
	    quicksort(a, index, left, i-1);
	    quicksort(a, index, i+1, right);
	}

	// partition a[left] to a[right], assumes left < right
	private static int partition(double[] a, double[] index, 
	int left, int right) {
	    int i = left - 1;
	    int j = right;
	    while (true) {
	        while (less(a[++i], a[right]))      // find item on left to swap
	            ;                               // a[right] acts as sentinel
	        while (less(a[right], a[--j]))      // find item on right to swap
	            if (j == left) break;           // don't go out-of-bounds
	        if (i >= j) break;                  // check if pointers cross
	        exch(a, index, i, j);               // swap two elements into place
	    }
	    exch(a, index, i, right);               // swap with partition element
	    return i;
	}

	// is x < y ?
	private static boolean less(double x, double y) {
	    return (x < y);
	}
	
	// exchange a[i] and a[j]
	private static void exch(double[] a, double[] index, int i, int j) {
	    double swap = a[i];
	    a[i] = a[j];
	    a[j] = swap;
	    double b = index[i];
	    index[i] = index[j];
	    index[j] = b;
	}
	
	public static double[][] computeBaseMass(String peptide, int[] pos, double[] ptmmass){
		double[] prefixMasses = new double[peptide.length()];
		double[] suffixMasses = new double [peptide.length()];
		double sum = 0.0; 
		//prefix sums
		int j = 0;
		for(int i = 0; i < peptide.length(); i++){
			if(pos != null && j < pos.length && pos[j]== i+1){
				sum += ptmmass[j];
				j++;
			}
			sum += Mass.getAAMass(peptide.charAt(i));
			prefixMasses[i] = sum;
			//System.out.println(sum);
		}
				
		j = 0;
		for(int i = 0; i < peptide.length()-1; i++){
			suffixMasses[i] = sum - prefixMasses[peptide.length()-i-2];
		}
		suffixMasses[peptide.length()-1] = 0; //we zero out the unfragmented precursor masses cause it 
		prefixMasses[peptide.length()-1] = 0; //do not provide information for IDs
		return new double[][]{prefixMasses, suffixMasses};		
	}
	/**Testing Code from here one===================================================================================================================================================**/
	public static void testGenerateTheoSpectrum(){
		String peptide = "ENEMLAQDK";
		String peptide2 = "ENEMLAQDK.2";
		long start = (new GregorianCalendar()).getTimeInMillis();
		Spectrum th = null;
		for(int i = 0; i < 100000; i ++){
			th = new TheoreticalSpectrum(peptide2);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + " created peaks: " + th.getPeak().size());
		start = (new GregorianCalendar()).getTimeInMillis();
		
		Spectrum s = null;
		for(int i = 0; i < 1000000; i ++){
			//s = getTheoSpectrum(peptide, 2, standardTypeMap, standardIonMap);
			s = getTheoSpectrumX(peptide, 2, standardTypeMap, standardIonMap);
		}
		System.out.println("matching " + 100 + " spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + " created peaks: " + s.getPeak().size());
	}
	
	public static void main(String[] args){
		testGenerateTheoSpectrum();
	}
	
}
