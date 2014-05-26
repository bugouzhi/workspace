//implement a special kind of spectrum that contain more than one spectrums
package org.Spectrums;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

public class MixtureSpectrum extends Spectrum{
	public int numberOfPeptide;
	public double[] spectrumCoefficient; //we view mixture as linear combination of single-peptide spectrum
	public Peptide[] peptides;
	
	public static List<Spectrum> createRandomMix(List<Spectrum> specList, int size, double minScale, double minSim, double maxSim, double massDiff, boolean toVector, int numSpect){
		int counts = 0;
		List<Spectrum> mixList = new ArrayList<Spectrum>();
		while(counts < size){
			mixList.add(createRandomMix(specList, minScale, minSim, maxSim, massDiff, toVector, numSpect));
			counts++;
		}
		return mixList;
	}
	
	public static Spectrum createRandomMix(List<Spectrum> specList, double minScale, double minSim, double maxSim, double massDiff, boolean toVector, int numSpect){
		List <Spectrum> mixtureSpect = new ArrayList<Spectrum>();; 
		Spectrum mixture, s1, s2;
		int counts = 0;
		int[] indexs = new int[numSpect];
		Spectrum[] selected = new Spectrum[numSpect];
		boolean pass = false;
		int numSelected = 0;
		//normalize all spectrum total intensity to one
		int count=0;
		while(numSelected < numSpect){
				int j = (int) (Math.random()*specList.size());
				s1 = specList.get(j);
				int k = 0;
				for(k = 0; k < numSelected; k++){
					s2 = selected[k];
					if((s1.cosineSim1(s2) < minSim) || (s1.cosineSim1(s2) > maxSim)  
							|| Math.abs(s1.parentMass - s2.parentMass)  > massDiff
							|| j == indexs[k]
							|| s1.charge > 3
							|| s1.peptide.contains("+")){
						count++;
						break;
					}
				}
				//System.out.println("k " + k);
				if(k==numSelected && s1.charge < 4 && !s1.peptide.contains("+")){
					indexs[numSelected] = j;
					selected[numSelected] = specList.get(j);
					numSelected++;
				}
				//if trying too many times we fail to find mixture we quit and retry
				if(count > 10000){
					return createRandomMix(specList, minScale, minSim, maxSim, massDiff, toVector, numSpect);
				}
		}
		for(int i = 0; i < selected.length; i++){
			selected[i].normalizeTIC();
		}
		Spectrum mix = selected[0];
		for(int i = 1; i < selected.length; i++){
			double scale = Math.random();
			scale = scale*(1.0-minScale)+minScale;
			selected[i].scaleSpectrum(scale);			
			mix = new Spectrum(mix, selected[i], toVector);	
		}
		mix.mergePeaks(mix, 0.3);
		//System.out.println("combined: " +  mix.getPeak().size());
		return mix;
	}
	
	public static void testCreatMix(){
		String specFile = "../mixture_linked/human_heck_1pepFDR_msgfdb.mgf";
		SpectrumLib lib = new SpectrumLib(specFile, "MGF");
		List<Spectrum> mix = createRandomMix(lib.getSpectrumList(), 1000, 0.3, 0.0, 0.4, 1.0, false, 3);
		SpectrumUtil.printSpectraToFile("../mixture_linked/3Mixtures.mgf", mix);
	}
	
	public static void main(String[] args){
		testCreatMix();
	}
	
	
	
	
    
}
