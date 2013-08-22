package org.Spectrums;

import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;
//contain  various method for analyzing mixture spectrum
//form long and compressed chromatography specifically
//for the Mixture Project

import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

public class ExperimentalMixtureAnalysis {
	
	//this method creat a graph the links distinct spectrums (as measure by overlap score) which
	//has similar mass and similar scan number, therefore are possible candidate to generate mixture
	//spectrum
	public static SimpleGraph getSpectrumMixtureGraph(SpectrumLib lib, double deltaM, double overlapCutOff, int deltaScans){
		Iterator<Spectrum> it = lib.iterator();
		Spectrum s = null;
		int counter = 0, i = 0, j=0;
		while(it.hasNext()){
			s = it.next();
			//System.out.println("name is: " + s.spectrumName);
			counter = Integer.parseInt((s.spectrumName.split("spec_|\\.")[1]));
			//System.out.println("counter is: " + counter);
			s.scanNumber = counter;
			s.score = s.scanNumber;
		}
		Collections.sort(lib.getSpectrumList());
		SimpleGraph<Spectrum, DefaultEdge> mixGraph = new SimpleGraph<Spectrum, DefaultEdge>(DefaultEdge.class);
		//every  spectrum is a vertex
		for(i = 0; i < lib.getSpectrumList().size(); i++){
			mixGraph.addVertex(lib.getSpectrumList().get(i));
		}
		
		Spectrum curr1, curr2, curr3;
		//an edge is created between two spectrum if they do not overlap, 
		//has small parentmass difference and has a scan number difference within
		//some range
		Set<Spectrum> vertices = mixGraph.vertexSet();
		Iterator<Spectrum> it1 = vertices.iterator();
		Iterator<Spectrum> it2;
		for(i = 0; i < lib.getSpectrumList().size(); i++){
			curr1 = lib.getSpectrumList().get(i);
			for(j = i+1; j < lib.getSpectrumList().size(); j++){
				curr2 = lib.getSpectrumList().get(j);
				if(Math.abs(curr1.scanNumber - curr2.scanNumber) < deltaScans){
					if(Math.abs(curr1.parentMass - curr2.parentMass) < deltaM){
						//System.out.println("overlap is: " + curr1.cosineSim(curr2));
						if(curr1.shiftCosineSim(curr2) < overlapCutOff){
							if(mixGraph.containsVertex(curr1) && mixGraph.containsVertex(curr2)){
								mixGraph.addEdge(curr1, curr2);    //need to make sure the vertex is not a duplicate and has not been remvoed
							}
						}else{
							//System.out.println(curr1.spectrumName + " overlap " + curr2.spectrumName + ": \t" + curr1.cosineSim(curr2));
							if(mixGraph.containsVertex(curr2)){
								if(curr1.sumOfPeaks() < curr2.sumOfPeaks()){
									mixGraph.removeVertex(curr2); //remove from the graph
								}else if(mixGraph.containsVertex(curr1)){
									mixGraph.removeVertex(curr1);
								}
							}
						}
					}
				}
			}
		}
		return mixGraph;
	}
	public static void testExperimentalData(){
		String filename = ".\\mixture_compressed\\new80min.mgf";
		String fileMix = ".\\mixture_compressed\\new2min.mgf";
		//String filename = ".\\MSPLib\\Lib\\human.msp";
        SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
		lib1.windowFilterPeaks(15, 50);
		mixlib.windowFilterPeaks(15, 50);
		lib1.printStat(SpectrumLib.NODETAIL);
		lib1.toNormVector(1, 0.5, 2000);
		mixlib.toNormVector(1, 0.5, 2000);
		
		//lib1.normIntensity();
		//mixlib.normIntensity();
		
		lib1.removeDuplicate(0.61);
		
		Spectrum temp = null;
		Vector v = null;
		Hashtable <String, Vector <Spectrum>> newTable = new Hashtable();
		int i;

		Iterator<Spectrum> it = mixlib.iterator();
		Spectrum mix, candidate, answer;
		Vector<Spectrum> spects;
		String[] peps, putativepeps;
		long start = (new GregorianCalendar()).getTimeInMillis();
		double accuracy = 0, size = mixlib.getAllSpectrums().size();
		i = 0;
		Vector<Spectrum> candidates;
		while(it.hasNext()){
			mix = it.next();
			if(mix.getPeak().size() == 0){ //make sure the mixture has some peaks otherwise it will take long time to search
				continue;
			}
			candidate = lib1.searchAndBoundLib(mix);
			//answer = lib1.topSpectrum(mix);
			//System.out.println(mix.peptide + ":\t" + answer.peptide + "\t" + answer.score + "\t" + mix.parentMass + "\t" + answer.parentMass);
		
		}
		System.out.println("Summary of Result: " + accuracy*100 + "%");
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
//	public static void statExperimentalData(){
//		String filename = ".\\mixture_compressed\\new80min.mgf";
//		String fileMix = ".\\mixture_compressed\\new2min.mgf";
//
//        SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
//		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
//		
//		lib1.windowFilterPeaks(15, 50);
//		mixlib.windowFilterPeaks(15, 50);
//		lib1.printStat(SpectrumLib.NODETAIL);
//		mixlib.printStat(SpectrumLib.NODETAIL);
//		lib1.toNormVector(1, 0.5, 2000);
//		mixlib.toNormVector(1, 0.5, 2000);
//		SimpleGraph sg = lib1.getNonRedundantSpectrumGraph(0.65);
//		lib1.getConnectedSpectrum(sg);
//		//SimpleGraph sgm = mixlib.getNonRedundantSpectrumGraph(0.55);
//		//lib1.getConnectedSpectrum(sgm);
//		
//	}
	
	public static void testSearchNBoundControl(){
		//String filename = "afumigatus_cmp_20.mgf";
		String filename = ".\\MSPLib\\Lib\\human.msp";
//		String filename = ".\\mixture_compressed\\new80min.mgf";
//		String fileMix = ".\\mixture_compressed\\new2min.mgf";

//        SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
        SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		lib1.removeModSpectra();
		lib1.printStat(SpectrumLib.NODETAIL);
		lib1.toNormVector(1, 0.5, 2000);
		//lib1.filterPeaks(20);
		SpectrumLib lib2 = lib1.Divide();
		Spectrum temp = null;
		Vector v = null;
		int i;
		SpectrumLib mixlib = lib2.createRandomMix(10000, 1, 0, 0.001, 1, 10000);
		Iterator<Spectrum> it = mixlib.iterator();
		Spectrum mix, candidate, answer;
		Vector<Spectrum> spects;
		String[] peps, putativepeps;
		long start = (new GregorianCalendar()).getTimeInMillis();
		double accuracy = 0, size = mixlib.getAllSpectrums().size();
		i = 0;
		List<Spectrum> spect1, spect2;
		while(it.hasNext()){
			mix = it.next();
			if(mix.getPeak().size() == 0){ //make sure the mixture has some peaks otherwise it will take long time to search
				continue;
			}
			System.out.println("searching matches to: " + mix.peptide);
			peps = mix.peptide.split(" & ");
			spect1 = lib1.getSpectra(peps[0]);
			spect2 = lib1.getSpectra(peps[1]);
//			lib1.spectrumLibrary.remove(peps[0]);
//			lib1.spectrumLibrary.remove(peps[1]);
//			lib1.spectrumList = lib1.getAllSpectrums();
			candidate = lib1.searchAndBoundLib(mix);
			//answer = lib1.findAnswerInSortedSpectrumList(mix, mix.peptide);
			//System.out.println(lib2.getSpectra(peps[0]).get(0));
			//System.out.println(lib2.getSpectra(peps[1]).get(0));
			putativepeps = candidate.peptide.split(" & ");
			if(peps[0].equals(putativepeps[0]) && peps[1].equals(putativepeps[1])){
				accuracy += 1/(size);
				System.out.println(candidate.peptide + " best score: " + candidate.score + "\t" + "t");
			}else if(peps[0].equals(putativepeps[1]) && peps[1].equals(putativepeps[0])){
				accuracy += 1/(size);
				System.out.println(candidate.peptide + " best score: " + candidate.score + "\t" + "t");
			}else{
				System.out.println(candidate.peptide + " best score: " + candidate.score + "\t" + "f");
			}
//			lib1.spectrumLibrary.put(peps[0], spect1);
//			lib1.spectrumLibrary.put(peps[1], spect2);
//			lib1.spectrumList = lib1.getAllSpectrums();
			i++;
			if(i == -1) {
				return;
			}
			
		}
		System.out.println("Summary of Result: " + accuracy*100 + "%");
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");

	}
	
	public static void testReverseSearchMix(){
		//String filename = ".\\mixture_compressed\\out1.mgf";
		//String fileMix = ".\\mixture_compressed\\out2.mgf";
		String filename = ".\\mixture_compressed\\new80min.mgf";		
		String fileMix = ".\\mixture_compressed\\new2min.mgf";
        //String fileMix = filename;
		//String fileMix = ".\\mixture_compressed\\newDirectInf.mgf";
        SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		lib1.printStat(SpectrumLib.NODETAIL);
		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
		mixlib.printStat(SpectrumLib.NODETAIL);
		
		lib1.windowFilterPeaks(15, 50);
		mixlib.windowFilterPeaks(15, 50);
		
		lib1.toNormVector(1, 0.5, 2000);
		mixlib.toNormVector(1, 0.5, 2000);
		
		lib1.normIntensity();
		mixlib.normIntensity();
		
		//get potential single-peptide spectrum that can form a mixture in
		//less-eluted samples
		SimpleGraph mixgraph = getSpectrumMixtureGraph(lib1, 2, 0.5, 1000);
		ConnectivityInspector conn = new ConnectivityInspector(mixgraph);
		List<Set<Spectrum>> c = conn.connectedSets();
		Iterator<Set<Spectrum>> comps = c.iterator(); //a list of connected components
		Object[] comp;
		//get all connected mixset has size 2 and 3 and try to search to see if they are presented
		//in the less eluded or direct infusion cases
		long start = (new GregorianCalendar()).getTimeInMillis();
		while(comps.hasNext()){
			comp = comps.next().toArray();
			if(comp.length <= 3){
				for(int i = 0; i < comp.length; i++){
					for(int j = i+1; j < comp.length; j++){
						mixlib.reverseSearchMix((Spectrum)comp[i], (Spectrum)comp[j]);
//						Spectrum s1 = (Spectrum)comp[i];
//						Spectrum s2 = (Spectrum)comp[j];
//						System.out.println(s1.peptide + " & " + s2.peptide + " with ratio: "
//								+ s1.sumOfPeaks() / s2.sumOfPeaks());
					}	
				}
			}
			
		}
		
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testReverseSearchAllPair(){
		String filename = ".\\mixture_compressed\\new80min.mgf";		
		String fileMix = ".\\mixture_compressed\\new2min.mgf";
        SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		lib1.printStat(SpectrumLib.NODETAIL);
		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
		mixlib.printStat(SpectrumLib.NODETAIL);
		Iterator<Spectrum> it = lib1.iterator();
		Spectrum s;
		int counter;
		while(it.hasNext()){
			s = it.next();
			//System.out.println("name is: " + s.spectrumName);
			counter = Integer.parseInt((s.spectrumName.split("spec_|\\.")[1]));
			//System.out.println("counter is: " + counter);
			s.scanNumber = counter;
			s.score = s.scanNumber;
		}
		lib1.windowFilterPeaks(15, 50);
		mixlib.windowFilterPeaks(15, 50);
		
		lib1.toNormVector(1, 0.5, 2000);
		mixlib.toNormVector(1, 0.5, 2000);
		
		lib1.normIntensity();
		mixlib.normIntensity();
		
		//get potential single-peptide spectrum that can form a mixture in
		//less-eluted samples
		//get all connected mixset has size 2 and 3 and try to search to see if they are presented
		//in the less eluded or direct infusion cases
		long start = (new GregorianCalendar()).getTimeInMillis();
		List<Spectrum> l = lib1.getAllSpectrums();
		Spectrum curr1, curr2;
		for(int i = 0; i < l.size(); i++ ){
			curr1 = l.get(i);
			for(int j = i + 1; j < l.size(); j++){
				curr2 = l.get(j);
				if(Math.abs(curr1.parentMass - curr2.parentMass) < 2){
					if(Math.abs(curr1.scanNumber - curr2.scanNumber) < 1000){
						if(curr1.cosineSim(curr2) < 0.5){
							mixlib.reverseSearchMix(curr1, curr2);
						}
					}
				}
			}
		}
		
		
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	//anohter control for revrerse search
	public static void testReverseSearchControl(){
		String filename = ".\\mixture_compressed\\out1.mgf";
		String fileMix = ".\\mixture_compressed\\out2.mgf";
        SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		lib1.printStat(SpectrumLib.NODETAIL);
		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
		mixlib.printStat(SpectrumLib.NODETAIL);

		//get potential single-peptide spectrum that can form a mixture in
		//less-eluted samples
		SimpleGraph mixgraph = getSpectrumMixtureGraph(lib1, 2, 0.5, 1000);
		ConnectivityInspector conn = new ConnectivityInspector(mixgraph);
		List<Set<Spectrum>> c = conn.connectedSets();
		Object[] comp1, comp2;
		//get all connected mixset has size 2 and 3 and try to search to see if they are presented
		//in the less eluded or direct infusion cases
		long start = (new GregorianCalendar()).getTimeInMillis();
		Spectrum s1, s2;
		for(int m = 0; m < c.size(); m++){
			for(int n = m+2; n < c.size(); n++){ //make sure they are two cluster away, so unlikely to be in same spectrum in compressed run
				comp1 = c.get(m).toArray();
				comp2 = c.get(n).toArray();
					for(int i = 0; i < 1; i++){
						for(int j = 0; j < 1; j++){
							s1 = (Spectrum)comp1[i];
							s2 = (Spectrum)comp2[j];
							if(Math.abs(s1.parentMass - s2.parentMass) <= 2){
								mixlib.reverseSearchMix((Spectrum)comp1[i], (Spectrum)comp2[j]);
							}
						}
					}
				}
			}
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	//control for reversesearch of mixture spectrum
	public static void testRandomRevSearchMix(){
		String filename = ".\\mixture_compressed\\out1.mgf";
		String fileMix = ".\\mixture_compressed\\out2.mgf";
		//String fileMix = ".\\mixture_compressed\\newDirectInf.mgf";
        
		SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		lib1.printStat(SpectrumLib.NODETAIL);
		//lib1.toNormVector(1, 0.5, 2000);
		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
		mixlib.printStat(SpectrumLib.NODETAIL);
		//mixlib.toNormVector(1, 0.5, 2000);
		int iterations = 1000;
		Spectrum s1, s2;
		long start = (new GregorianCalendar()).getTimeInMillis();
		while(iterations > 0){
			s1 = lib1.getRandomSpectrum();
			s2 = lib1.getRandomSpectrum();
			//System.out.println("overlap is: " + s1.shiftCosineSim(s2));
			//System.out.println("delta m is: " + Math.abs(s1.parentMass - s2.parentMass));
			if(s1.shiftCosineSim(s2) < 0.5 && Math.abs(s1.parentMass - s2.parentMass) > 5){
				mixlib.reverseSearchMix(s1, s2);
				iterations--;
			}
		}
		
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
}
