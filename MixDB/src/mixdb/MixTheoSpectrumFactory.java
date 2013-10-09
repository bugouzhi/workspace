package mixdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.Spectrums.Mass;
import org.Spectrums.Peptide;
import org.Spectrums.Spectrum;
import org.Spectrums.TheoreticalSpectrum;

/**
 * Theoretical spectrum factory for mixture spectrum
 * @author Jian Wang
 *
 */
public class MixTheoSpectrumFactory{
	private static int minCharge = 2;
	private static int maxCharge = 3;
	public static int MixTypeOffSet=0;
	public static Map<String, MixIonType> mixTypeMap=createStandardMixMap();
	public static IonTypeMapper mixIonMap;
	public static Map<String, Collection<IonType>> mixPeptideMap;
	
	public static Map<String, MixIonType> createStandardMixMap(){
		Map<String, IonType> standardMap = TheoreticalSpectrumFactory.standardTypeMap;
		System.out.println("standard ions counts:  " + standardMap.values().size());
		Map<String, MixIonType> mixTypeMap = new HashMap<String, MixIonType>();
		List<IonType> mix = new ArrayList<IonType>();
		List<IonType> mixLow = new ArrayList<IonType>();
		int count=0;
		for(Iterator<String> it = standardMap.keySet().iterator(); it.hasNext();){
			String key = it.next();
			for(int c = minCharge; c <= maxCharge; c++){
				SimpleIonType type = (SimpleIonType)standardMap.get(key);
				int mixtureCharge = type.getPepCharge()+c;
				MixIonType hi = new MixIonType(type, mixtureCharge, MixturePeptideType.HIGH_ABUNDANCE);
				MixIonType low = new MixIonType(type, mixtureCharge, MixturePeptideType.LOW_ABUNDANCE);
				mixTypeMap.put(key+"@"+mixtureCharge+"@high", hi);
				mixTypeMap.put(key+"@"+mixtureCharge+"@low", low);
				mix.add(hi);
				mixLow.add(low);
				count+=2;
			}
		}
		MixTheoSpectrumFactory.MixTypeOffSet = mix.size();
		mix.addAll(mixLow);
		MixTheoSpectrumFactory.mixIonMap = IonTypeMapper.createIonTypeMap(mix);
		MixTheoSpectrumFactory.mixPeptideMap = createMixPeptideMap(mix);
		MixTheoSpectrumFactory.mixTypeMap = mixTypeMap;
		System.out.println("standard mix ions counts:  " + mixTypeMap.size());
		System.out.println("standard mix ions counts:  " + count);
		return mixTypeMap;
	}
	
	/**
	 * transform a set of simple ion type to mix ion types
	 * @param type
	 * @return
	 */
	public static List<MixIonType> getMixIonType(List<SimpleIonType> types){
		List<MixIonType> mix = new ArrayList<MixIonType>();
		List<MixIonType> mixLow = new ArrayList<MixIonType>();
		for(Iterator<SimpleIonType> it = types.iterator(); it.hasNext();){
			for(int c = minCharge; c <= maxCharge; c++){
				SimpleIonType type = it.next();
				int mixtureCharge = type.getPepCharge()+c;
				MixIonType hi = new MixIonType(type, mixtureCharge, MixturePeptideType.HIGH_ABUNDANCE);
				MixIonType low = new MixIonType(type, mixtureCharge, MixturePeptideType.LOW_ABUNDANCE);
				mix.add(hi);
				mixLow.add(low);
			}
		}
		mix.addAll(mixLow);
		return mix;
	}
	
	public static Map<String, Collection<IonType>> createMixPeptideMap(List<IonType> ionTypeList){
		Map<String, Collection<IonType>> peptideMap = new HashMap<String, Collection<IonType>>();
		for(Iterator<IonType> it = ionTypeList.iterator(); it.hasNext();){
			MixIonType type = (MixIonType)it.next();
			MixturePeptideType pType = type.getPeptideType();
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
			if(p.length() <= 18){
				return SimplePeptideType.SHORT;
			}else{
				return SimplePeptideType.LONG;
			}
		}
		return 0;
	}
	
	public static Spectrum getMixTheoSpectrum(ArrayTheoreticalSpectrum th1, ArrayTheoreticalSpectrum th2){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] massIntList1 = transformIonType(th1.getMassIntensityList(), th2.getCharge(), 0);//MixturePeptideType.HIGH_ABUNDANCE);
		double[][] massIntList2 = transformIonType(th2.getMassIntensityList(), th1.getCharge(), 1);//MixturePeptideType.LOW_ABUNDANCE);
		double[][] massIntList = merge(massIntList1, massIntList2);
		//System.out.println("size: " + massIntList[0].length);
		arry.setMassIntensityList(massIntList);
		arry.peptide = th1.peptide + " & " + th2.peptide;
		arry.parentMass = th1.parentMass;
		arry.charge = th1.charge + th2.charge;
		return arry;
	}
	
	//now we add an nonmix to mix
	public static Spectrum getMixTheoSpectrum(ArrayTheoreticalSpectrum mix, ArrayTheoreticalSpectrum th2, int mixNum){
		if(mixNum == 2){
			return getMixTheoSpectrum(mix, th2);
		}
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] massIntList1 = mix.getMassIntensityList();//MixturePeptideType.HIGH_ABUNDANCE);
		double[][] massIntList2 = transformIonType(th2.getMassIntensityList(), 2, 1);//we can no longer do mix charge here since it is unknown
		double[][] massIntList = merge(massIntList1, massIntList2);
		//System.out.println("size: " + massIntList[0].length);
		arry.setMassIntensityList(massIntList);
		arry.peptide = mix.peptide + " & " + th2.peptide;
		arry.parentMass = mix.parentMass;
		arry.charge = mix.charge;
		return arry;
	}
	
	
	public static double[][] transformIonType(double[][] massIntList, int partnerCharge, int abundance){
		double[][] transformed = new double[2][massIntList[0].length];
		for(int i = 0; i < massIntList[0].length; i++){
			transformed[0][i] = massIntList[0][i]; 
			transformed[1][i] = getMixIonType(massIntList[1][i], partnerCharge, abundance);
		}
		return transformed;
	}
	
	public static double getMixIonType(double type, int partnerCharge, int abundance){
		//System.out.println("original type: " + TheoreticalSpectrumFactory.standardIonMap.getIonType((int)type) +"\t" + type);
		int index = MixTheoSpectrumFactory.maxCharge-partnerCharge;
		double offset = type*(maxCharge-minCharge+1) + (partnerCharge-maxCharge)+2;
		double mixtype = offset-1+abundance*MixTheoSpectrumFactory.MixTypeOffSet;
		//System.out.println("mix type: " + mixIonMap.getIonType((int)mixtype) + "\t" + mixtype);
		return mixtype;
	}
	
	public static Spectrum getMixTheoSpectrum(String peptide1, String peptide2, int charge1, int charge2, Map<String, MixIonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] massIntList1 = getMixMassIntList(peptide1, charge1, charge1+charge2, "high", typeMap, ionMap);
		double[][] massIntList2 = getMixMassIntList(peptide2, charge2, charge1+charge2, "low", typeMap, ionMap);
		double[][] massIntList = merge(massIntList1, massIntList2);
		Peptide p1 = new Peptide(peptide1+"."+charge1);
		Peptide p2 = new Peptide(peptide2+"."+charge2);
		//System.out.println("size: " + massIntList[0].length);
		arry.setMassIntensityList(massIntList);
		arry.peptide = peptide1+"." + charge1 + " & " + peptide2 + "." + charge2;
		arry.parentMass = p1.getParentmass();
		arry.charge = charge1+charge2;
		return arry;
	}
	
	public static double[][] merge(double[][] massIntList1, double[][] massIntList2){
		double[][] massIntList = new double[2][massIntList1[0].length+massIntList2[0].length];
		int i = 0, j = 0;
		int curr = 0;
		while(i < massIntList1[0].length && j < massIntList2[0].length){
			if(massIntList1[0][i] <= massIntList2[0][j]){
				massIntList[0][curr] = massIntList1[0][i];
				massIntList[1][curr] = massIntList1[1][i];
				curr++;
				i++;
			}else{
				massIntList[0][curr] = massIntList2[0][j];
				massIntList[1][curr] = massIntList2[1][j];
				curr++;
				j++;
			}
		}
		while(i < massIntList1[0].length){
			massIntList[0][curr] = massIntList1[0][i];
			massIntList[1][curr] = massIntList1[1][i];
			curr++;
			i++;
		}
		while(j < massIntList2[0].length){
			massIntList[0][curr] = massIntList2[0][j];
			massIntList[1][curr] = massIntList2[1][j];
			curr++;
			j++;
		}
		return massIntList;
	}
	
	public static double[][] getMixMassIntList(String peptide, int charge, int mixtureCharge, String abundance, Map<String, MixIonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] baseMasses = TheoreticalSpectrumFactory.computeBaseMass(peptide, new int[]{}, new double[]{}); 
		String[] leng = new String[]{"short", "long"};
		Collection<IonType> types = mixPeptideMap.get(""+charge+"@"+leng[getPeptideLength(peptide, charge)]+"@"+mixtureCharge+"@"+abundance);
		double[][] massIntList = new double[2][baseMasses[0].length*types.size()];
		int currentIndex = 0;
		for(Iterator<IonType> it = types.iterator(); it.hasNext();){
			SimpleIonType type = (SimpleIonType)it.next();
			int index = mixIonMap.getIndex(type);
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
		TheoreticalSpectrumFactory.quicksort(massIntList[0], massIntList[1]);
		return massIntList;
	}
	
}
