package mixdb;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.Spectrums.LinkedPeptide;
import org.Spectrums.Mass;
import org.Spectrums.Peptide;
import org.Spectrums.Spectrum;
import org.Spectrums.TheoreticalSpectrum;

public class LinkedTheoSpectrumFactory {
	private static int minCharge = 0;
	private static int maxCharge = 3;
	public static int LinkedTypeOffSet=0;
	public static int MixOffSet = 0;
	public static Map<String, LinkedIonType> linkedTypeMap=createStandardLinkMap();
	public static IonTypeMapper linkedIonMap;
	public static Map<String, Collection<IonType>> linkedPeptideMap;
	
	public static Map<String, LinkedIonType> createStandardLinkMap(){
		String[] pTypes = TheoreticalSpectrum.prefixIons;
		String[] sTypes = TheoreticalSpectrum.suffixIons;
		int[] peptideCharge = {2,3,4,5};
		int[] peptideLength = {SimplePeptideType.SHORT}; //no length modelling for linked peptides here
		Map<String, IonType> standardMap = TheoreticalSpectrumFactory.createIonTypeMap(pTypes, sTypes, peptideCharge, peptideLength);
		System.out.println("standard ions counts:  " + standardMap.values().size());
		Map<String, LinkedIonType>  linkedTypeMap = new HashMap<String, LinkedIonType>();
		List<IonType> linkedHi = new ArrayList<IonType>();
		List<IonType> unLinkedHi = new ArrayList<IonType>();
		List<IonType> linkedLo = new ArrayList<IonType>();
		List<IonType> unLinkedLo = new ArrayList<IonType>();
		List<IonType> linked = new ArrayList<IonType>();
		int count=0;
		//creating linked ion type
		for(Iterator<String> it = standardMap.keySet().iterator(); it.hasNext();){
			String key = it.next();
			SimpleIonType type = (SimpleIonType)standardMap.get(key);
			LinkedIonType linkTypeHi = new LinkedIonType(new MixIonType(type, type.getPepCharge(), MixturePeptideType.HIGH_ABUNDANCE), true);
			LinkedIonType unlinkTypeHi = new LinkedIonType(new MixIonType(type, type.getPepCharge(), MixturePeptideType.HIGH_ABUNDANCE), false);
			LinkedIonType linkTypeLo = new LinkedIonType(new MixIonType(type, type.getPepCharge(), MixturePeptideType.LOW_ABUNDANCE), true);
			LinkedIonType unlinkTypeLo = new LinkedIonType(new MixIonType(type, type.getPepCharge(), MixturePeptideType.LOW_ABUNDANCE), false);
			linkedTypeMap.put(key+"@Hi"+"@linked", linkTypeHi);
			linkedTypeMap.put(key+"@Hi"+"@unlinked", unlinkTypeHi);
			linkedTypeMap.put(key+"@Lo"+"@linked", linkTypeLo);
			linkedTypeMap.put(key+"@Lo"+"@unlinked", unlinkTypeLo);
			linkedHi.add(linkTypeHi);
			linkedLo.add(linkTypeLo);
			count+=2;
			unLinkedHi.add(unlinkTypeHi);
			unLinkedLo.add(unlinkTypeLo);
			count+=2;
		}
		LinkedTheoSpectrumFactory.LinkedTypeOffSet = linkedHi.size();
		LinkedTheoSpectrumFactory.MixOffSet = unLinkedHi.size() + linkedHi.size();
		linked.addAll(unLinkedHi);
		linked.addAll(linkedHi);
		linked.addAll(unLinkedLo);
		linked.addAll(linkedLo);
		LinkedTheoSpectrumFactory.linkedIonMap = IonTypeMapper.createIonTypeMap(linked);
		LinkedTheoSpectrumFactory.linkedPeptideMap = createLinkedPeptideMap(linked);
		//LinkedTheoSpectrumFactory.linkedTypeMap = linkedTypeMap;
		System.out.println("standard linked ions counts:  " + linked.size());
		System.out.println("standard linked ions counts:  " + count);
		return linkedTypeMap;
	}
	
	
	public static Map<String, Collection<IonType>> createLinkedPeptideMap(List<IonType> ionTypeList){
		Map<String, Collection<IonType>> peptideMap = new HashMap<String, Collection<IonType>>();
		for(Iterator<IonType> it = ionTypeList.iterator(); it.hasNext();){
			LinkedIonType type = (LinkedIonType)it.next();
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
				if(!type.isLinked())         //only added unlinked type... we shift to linked type when see a linked fragment
					types.add(type);
			}
			peptideMap.put(key, types);
		}
		for(Iterator<String> it = peptideMap.keySet().iterator(); it.hasNext();){
			String key = it.next();
			Collection<IonType> types = peptideMap.get(key);
			TreeMap<Double, IonType> sortedType = new TreeMap();
			//System.out.println("we start with types: " + types.size());
			for(Iterator<IonType> it2 = types.iterator(); it2.hasNext();){
				SimpleIonType curr = (SimpleIonType)it2.next();
				double mass = (curr.getOffset()+500)/curr.getCharge() + 5000*curr.getDirection(); //try to separate ion type
				//System.out.println(mass);
				sortedType.put(mass, curr);
			}
			types.clear();
			while(!sortedType.isEmpty()){
				Map.Entry<Double, IonType> entry = sortedType.pollFirstEntry();
				//System.out.println("type: " + entry.getKey() + "\t" + entry.getValue());
				types.add(entry.getValue());
			}
			//System.out.println("sorting and with types: " + types.size());
			

		}
		return peptideMap;
	}
	
	public static int getPeptideLength(Peptide p){
		return getPeptideLength(p.getPeptide(), p.getCharge());
	}
	
	public static int getPeptideLength(String p, int charge){
//		if(charge == 2){
//			if(p.length() <= 12){
//				return SimplePeptideType.SHORT;
//			}else{
//				return SimplePeptideType.LONG;
//			}
//		}else if(charge == 3){
//			if(p.length() <= 18){
//				return SimplePeptideType.SHORT;
//			}else{
//				return SimplePeptideType.LONG;
//			}
//		}
		//noted: there is no model of length in linked peptides
		return SimplePeptideType.SHORT;
	}
	
	public static Spectrum getMixTheoSpectrum(ArrayTheoreticalSpectrum th1, ArrayTheoreticalSpectrum th2){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] massIntList1 = transformIonType(th1.getMassIntensityList(), th2.getCharge(), 0);//MixturePeptideType.HIGH_ABUNDANCE);
		double[][] massIntList2 = transformIonType(th2.getMassIntensityList(), th1.getCharge(), 1);//MixturePeptideType.LOW_ABUNDANCE);
		double[][] massIntList = MixTheoSpectrumFactory.merge(massIntList1, massIntList2);
		//System.out.println("size: " + massIntList[0].length);
		arry.setMassIntensityList(massIntList);
		arry.peptide = th1.peptide + " & " + th2.peptide;
		arry.parentMass = th1.parentMass;
		arry.charge = th1.charge + th2.charge;
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
		int index = LinkedTheoSpectrumFactory.maxCharge-partnerCharge;
		double offset = type*(maxCharge-minCharge+1) + (partnerCharge-maxCharge)+2;
		double mixtype = offset-1+abundance*MixTheoSpectrumFactory.MixTypeOffSet;
		//System.out.println("mix type: " + mixIonMap.getIonType((int)mixtype) + "\t" + mixtype);
		return mixtype;
	}
	
	
	public static Spectrum getLinkedTheoSpectrum(String peptide, int charge, int[] pos, double[] ptmmasses,  int linkedPos, 
			Map<String, LinkedIonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] massIntList = getLinkedMassIntList(peptide, charge, charge, pos, ptmmasses, linkedPos, "high", typeMap, ionMap);
		Peptide p1 = new Peptide(peptide+"."+charge);
		//System.out.println("size: " + massIntList[0].length);
		arry.setMassIntensityList(massIntList);
		arry.peptide = peptide+"." + charge;
		arry.parentMass = p1.getParentmass();
		arry.charge = charge;
		return arry;
	}
	
	
	public static Spectrum getLinkedTheoSpectrum(LinkedPeptide lp, Map<String, LinkedIonType> typeMap, IonTypeMapper ionMap){
		Peptide[] peps = lp.getPeptides();
		Peptide p1 = peps[0];
		Peptide p2 = peps[1];
		return getLinkedTheoSpectrum(p1.getPeptide(), p2.getPeptide(), p1.getCharge(), 
				p2.getCharge(), p1.getPos(), p2.getPos(), p1.getPtmmasses(), p2.getPtmmasses(),
				p1.getLinkedPos(), p2.getLinkedPos(), typeMap, ionMap);
	}
	
	public static Spectrum getLinkedTheoSpectrum(String peptide1, String peptide2, int charge1, int charge2, 
			int[] pos1, int[] pos2, double[] ptmmasses1, double[] ptmmasses2, int linkedPos1, int linkedPos2,
			Map<String, LinkedIonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] massIntList1 = getLinkedMassIntList(peptide1, charge1, charge1, pos1, ptmmasses1, linkedPos1, "high", typeMap, ionMap);
		double[][] massIntList2 = getLinkedMassIntList(peptide2, charge2, charge2, pos2, ptmmasses2, linkedPos2, "low", typeMap, ionMap);
		double[][] massIntList = MixTheoSpectrumFactory.merge(massIntList1, massIntList2);
		Peptide p1 = new Peptide(peptide1+"."+charge1);
		Peptide p2 = new Peptide(peptide2+"."+charge2);
		//System.out.println("size: " + massIntList[0].length);
		arry.setMassIntensityList(massIntList);
		arry.peptide = peptide1+"." + charge1 + " & " + peptide2 + "." + charge2;
		arry.parentMass = p1.getParentmass();
		arry.charge = charge1+charge2;
		return arry;
	}
	
	public static int getLinkedIndex(MixIonType type){
		return linkedIonMap.getIndex(type)+LinkedTheoSpectrumFactory.LinkedTypeOffSet;
	}
	
	public static double[][] getLinkedMassIntList(String peptide, int charge, int mixtureCharge, 
			int[] pos, double[] ptmmasses, int linkedPos, String abundance, 
			Map<String, LinkedIonType> typeMap, IonTypeMapper ionMap){
		if(charge > 5){
			charge = 5;
		}
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] baseMasses = TheoreticalSpectrumFactory.computeBaseMass(peptide, pos, ptmmasses); 
		String[] leng = new String[]{"short", "long"};
		Collection<IonType> types = linkedPeptideMap.get(""+charge+"@"+leng[getPeptideLength(peptide, charge)]+"@"+charge+"@"+abundance);
		double[][] massIntList = new double[2][baseMasses[0].length*types.size()];
		int currentIndex = 0;
//		System.out.println(pos[pos.length-1] + "\t" + ptmmasses[pos.length-1] +"\t" + linkedPos);
		for(Iterator<IonType> it = types.iterator(); it.hasNext();){
			SimpleIonType type = (SimpleIonType)it.next();
			int index = linkedIonMap.getIndex(type);
			int linkedIndex = getLinkedIndex((MixIonType)type);
//			if(linkedIndex > linkedIonMap.getSize()){
//				System.out.println("ion type out of range");
//			}
			if(type.getDirection() == SimpleIonType.PREFIX){
				for(int i = 0; i < linkedPos-1; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntList[0][currentIndex] = mass;
					massIntList[1][currentIndex++] = index;
				}
				//
				for(int i = linkedPos-1; i < baseMasses[0].length; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntList[0][currentIndex] = mass;
					massIntList[1][currentIndex++] = linkedIndex;
				}
			}		
			if(type.getDirection() == SimpleIonType.SUFFIX){
				int link = baseMasses[1].length-linkedPos;
				for(int i = 0; i < link; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[1][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntList[0][currentIndex] = mass;
					massIntList[1][currentIndex++] = index;
				}
				
				for(int i = link; i < baseMasses[1].length; i++){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[1][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntList[0][currentIndex] = mass;
					massIntList[1][currentIndex++] = linkedIndex;
				}
			}
		}
		TheoreticalSpectrumFactory.quicksort(massIntList[0], massIntList[1]);
		return massIntList;
	}
	
	//a faster version
	public static double[][] getLinkedMassIntList2(String peptide, int charge, int mixtureCharge, 
			int[] pos, double[] ptmmasses, int linkedPos, String abundance, 
			Map<String, LinkedIonType> typeMap, IonTypeMapper ionMap){
		ArrayTheoreticalSpectrum arry = new ArrayTheoreticalSpectrum();
		double[][] baseMasses = TheoreticalSpectrumFactory.computeBaseMass(peptide, pos, ptmmasses); 
		String[] leng = new String[]{"short", "long"};
		Collection<IonType> types = linkedPeptideMap.get(""+charge+"@"+leng[getPeptideLength(peptide, charge)]+"@"+mixtureCharge+"@"+abundance);
		double[][] massIntList = new double[2][baseMasses[0].length*types.size()];
		double[][] massIntListPrefix = new double[2][baseMasses[0].length*types.size()/2];
		double[][] massIntListSuffix = new double[2][baseMasses[0].length*types.size()/2];
		int currentIndex = 0;
//		System.out.println(pos[pos.length-1] + "\t" + ptmmasses[pos.length-1] +"\t" + linkedPos);
		
		//first do prefix
		for(int i = 0; i < linkedPos-1; i++){
			for(Iterator<IonType> it = types.iterator(); it.hasNext();){
				SimpleIonType type = (SimpleIonType)it.next();
				int index = linkedIonMap.getIndex(type);
				if(type.getDirection() == SimpleIonType.PREFIX){
					//System.out.println("type" + type.toString());
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					//System.out.println("mass " + mass);
					massIntListPrefix[0][currentIndex] = mass;
					massIntListPrefix[1][currentIndex++] = index;
				}
			}
		}
		System.out.println("prefix done");
		for(int i = linkedPos-1; i < baseMasses[0].length; i++){
			for(Iterator<IonType> it = types.iterator(); it.hasNext();){
				SimpleIonType type = (SimpleIonType)it.next();
				int linkedIndex = getLinkedIndex((MixIonType)type);
				if(type.getDirection() == SimpleIonType.PREFIX){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntListPrefix[0][currentIndex] = mass;
					massIntListPrefix[1][currentIndex++] = linkedIndex;
				}
			}
		}
		System.out.println("suffix done");
		//then do suffix
		currentIndex = 0; //reset the index
		int link = baseMasses[1].length-linkedPos;
		for(int i = 0; i < link; i++){
			for(Iterator<IonType> it = types.iterator(); it.hasNext();){
				SimpleIonType type = (SimpleIonType)it.next();
				int index = linkedIonMap.getIndex(type);
				if(type.getDirection() == SimpleIonType.SUFFIX){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntListSuffix[0][currentIndex] = mass;
					massIntListSuffix[1][currentIndex++] = index;
				}
			}
		}
		
		for(int i = link; i < baseMasses[1].length; i++){
			for(Iterator<IonType> it = types.iterator(); it.hasNext();){
				SimpleIonType type = (SimpleIonType)it.next();
				int linkedIndex = getLinkedIndex((MixIonType)type);
				if(type.getDirection() == SimpleIonType.SUFFIX){
					int peakCharge = type.getCharge();
					double mass = (baseMasses[0][i]+type.getOffset()+Mass.PROTON_MASS*(peakCharge-1))/peakCharge;
					massIntListSuffix[0][currentIndex] = mass;
					massIntListSuffix[1][currentIndex++] = linkedIndex;
				}
			}
		}
		int i = 0, j = 0;
		currentIndex = 0;
		while(i < massIntListPrefix[0].length && j < massIntListSuffix[0].length){
			if(massIntListPrefix[0][i] < massIntListSuffix[0][j]){
				massIntList[0][currentIndex] = massIntListPrefix[0][i];
				massIntList[1][currentIndex] = massIntListPrefix[1][i];
				i++;
			}else{
				massIntList[0][currentIndex] = massIntListSuffix[0][j];
				massIntList[1][currentIndex] = massIntListSuffix[1][j];
				j++;
			}
			currentIndex++;
		}
		for(;i < massIntListPrefix[0].length; i++){
			massIntList[0][currentIndex] = massIntListPrefix[0][i];
			massIntList[1][currentIndex] = massIntListPrefix[1][i];
			currentIndex++;
		}
		for(;j < massIntListSuffix[0].length; j++){
			massIntList[0][currentIndex] = massIntListSuffix[0][j];
			massIntList[1][currentIndex] = massIntListSuffix[1][j];
			currentIndex++;
		}
		
		//TheoreticalSpectrumFactory.quicksort(massIntList[0], massIntList[1]);
		return massIntList;
	}

}
