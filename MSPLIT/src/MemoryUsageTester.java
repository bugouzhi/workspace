

import java.util.ArrayList;
import java.util.List;

import Spectrum.Peak;
import Spectrum.SpectrumLib;

public class MemoryUsageTester {
	public static void testLoadSpectrum(){
		//String file = "..\\MSPLib\\Lib\\yeast.msp";
		String file = "../mixture_linked/spectral_library/NIST_yeast_IT_v3.0_2009-05-04_7AA_decoy.sptxt";
		getMemoryUsage();
		//SpectrumLib lib1 = new SpectrumLib(file, "MSP");
		SpectrumLib lib1 = new SpectrumLib(file, "splib");
		getMemoryUsage();
		Runtime.getRuntime().gc();
		getMemoryUsage();
		System.out.println(lib1.getAllSpectrums().size());
	}
	
	public static void testLoadPeaks(){
		getMemoryUsage();
		List<Peak> pList = new ArrayList<Peak>();
		for(int i = 0; i < 5425265; i++){
			Peak p = new Peak(Math.random(), Math.random());
			pList.add(p);
		}
		getMemoryUsage();
		Runtime.getRuntime().gc();
		getMemoryUsage();
		System.out.println(pList.size());
	}
	
	public static void getMemoryUsage(){
		Runtime r = Runtime.getRuntime();
		System.out.println("Current memory usage: " 
				+ (r.totalMemory()
						- r.freeMemory()));
	}
	
	public static void main(String[] args){
		//testLoadPeaks();
		testLoadSpectrum();
	}
}
