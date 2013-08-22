package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

public class MemoryUsageTester {
	public static void testLoadPeptides(){
		String file = "..\\mixture_linked\\database\\Yeast_allPeptides_plusDecoy.txt";
		getMemoryUsage();
		List<String> l = Utils.FileIOUtils.createListFromFile(file);
		getMemoryUsage();
		List list = new ArrayList();
		StringBuffer protein = new StringBuffer();
		int begin=0;
//		for(int i = 0; i < 1500000; i++){
//			protein.append(l.get(i));
//			list.add(new PeptideLite(begin, protein.length(), "", 1));
//		}
		
		for(int i = 0; i < 1500000; i++){
			//Peptide p1 = new Peptide(l.get(i), 1);
			//list.add(p1);
			//list.add(p1);
			//list.add(p2);
			//list.add(new Peptide("",2));
			//list.add(new Peptide("",3));
		}
		
		for(int i = 0; i < 854; i ++){
			Peptide p1 = new Peptide(l.get(i), 1);
			for(int j = 0; j < 855; j++){
				Peptide p2 = new Peptide(l.get(j), 1);
				for(int c = 3; c < 7; c++){
					LinkedPeptide lp = new LinkedPeptide(p1, p2, c, 1, 1);
					list.add(lp);
				}
			}
		}
		getMemoryUsage();
		Runtime.getRuntime().gc();
		getMemoryUsage();
		System.out.println(list.size());
		//System.out.println(longStr.length());
	}
	
	public static void testLoadSpectrum(){
		String file = "..\\MSPLib\\Lib\\ecoli.msp";
		getMemoryUsage();
		SpectrumLib lib1 = new SpectrumLib(file, "MSP");
		getMemoryUsage();
		Runtime.getRuntime().gc();
		getMemoryUsage();
		System.out.println(lib1.getAllSpectrums().size());
	}
	
	public static void getMemoryUsage(){
		Runtime r = Runtime.getRuntime();
		System.out.println("Current memory usage: " 
				+ (r.totalMemory()
						- r.freeMemory()));
	}
	
	public static void main(String[] args){
		testLoadPeptides();
	}
}
