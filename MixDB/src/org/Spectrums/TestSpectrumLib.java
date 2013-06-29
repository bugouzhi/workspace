package org.Spectrums;

import java.util.GregorianCalendar;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

/**
 * Various method for testing SpectrumLib 
 * @author jian wang
 *
 */
public class TestSpectrumLib {
	public static void testWeightedCosineSim(){
		String filename =".\\MSPLib\\Lib\\human.msp";
		SpectrumLib l = new SpectrumLib(filename, "MSP") ;
		Spectrum mix ;
        
		Iterator it = l.getSpectrumLibrary().values().iterator();
		
		l.removeModSpectra() ;
		l.removeSingle() ;
		l.toNormVector(1, 0.5, 2000);
		SpectrumLib l2 = l.Divide() ;
		SpectrumLib lib1 = l2.createMix("bin2.name", 1000, 1, 0.2, 0, 1, 2000) ;
		System.out.println("after dividing into two");
		lib1.toNormVector(1, 0.5, 2000);
		l.normIntensity();
		lib1.normIntensity();
		//l.toNormVector(1, 0.5, 2000);
		//lib1.toNormVector(1, 0.5, 2000);
		it = l.iterator();	
		List temp = lib1.getAllSpectrums() ;
		String peptides[];
		Spectrum a, b ;
		//positive sets
		for (int j = 0; j < 999; j++) {
			double alpha;
			mix = (Spectrum)temp.get(j) ;
			peptides = mix.peptide.split(" & ") ;
			a = l.getSpectra(peptides[0]).get(0) ;
			b = l.getSpectra(peptides[1]).get(0) ;
			alpha = mix.alpha(a, b);
			System.out.println(mix.peptide + "\tAlpha: " + mix.alpha(a, b) + "\tScore:" + mix.maxScore(a, b, alpha) + "\t" + mix.maxScore(a,b, 2)) ;
			//System.out.println("Alpha: " + mix.alpha(a, b) + "\tScore:" + mix.maxScore(a, b, 0.5)) ;
			
		}
	
	}
	
	public static void testResidualAlpha(){ //test estimation of alpha by using the residual intensity of a mxiture
		String filename =".\\MSPLib\\Lib\\human.msp";
		SpectrumLib l = new SpectrumLib(filename, "MSP") ;
		Spectrum mix ;
		Iterator it = l.getSpectrumLibrary().values().iterator();
		l.removeModSpectra() ;
		l.removeSingle() ;
		//l.windowFilterPeaks(5, 50);
		l.toNormVector(1, 0.5, 2000);
		SpectrumLib l2 = l.Divide() ;
		SpectrumLib lib1 = l2.createRandomMix(1000, 1, 1, 0, 1, 2000) ;
		lib1.toNormVector(1, 0.5, 2000);
		lib1.normIntensity();
		l.normIntensity();
		it = l.iterator();	
		List temp = lib1.getAllSpectrums() ;
		String peptides[];
		Spectrum a, b ;
		//positive sets
		for (int j = 0; j < 1000; j++) {
			double alpha;
			mix = (Spectrum)temp.get(j) ;
			peptides = mix.peptide.split(" & ") ;
			a = l.getSpectra(peptides[0]).get(0) ;
			b = l.getSpectra(peptides[1]).get(0) ;
			alpha = mix.residual(a);
			//alpha = mix.residualIntensity(a);
			System.out.println("Alpha: " + mix.residual(a) + "\tScore:" + mix.maxScore(a, b, alpha) + "\t" + mix.maxScore(a,b, 0.1)) ;
			//System.out.println("Alpha: " + mix.alpha(a, b) + "\tScore:" + mix.maxScore(a, b, 0.5)) ;
			
		}
	
	}
	
	public static void testReadMGF(){
		String filename = "afumigatus_cmp_20.mgf";
		SpectrumLib l = new SpectrumLib(filename, "MGF");
		System.out.print(l);
		l.printStat(SpectrumLib.DETAIL);
		l.removeModSpectra();
		l.printStat(SpectrumLib.NODETAIL);
		System.out.println();
		System.out.println("Dividing Spectrums library into two sets");
		SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		lib1.printStat(SpectrumLib.NODETAIL);
		lib1.removeModSpectra();
		lib1.printStat(SpectrumLib.NODETAIL);
		SpectrumLib lib2 = lib1.Divide();
		lib1.printStat(SpectrumLib.NODETAIL);
		lib2.printStat(SpectrumLib.NODETAIL);
		SpectrumLib mixlib = lib2.createMix(5);
		mixlib.printStat(SpectrumLib.NODETAIL);
		l.getSpectrumSimilarity(l.getSpectrumList().size());
		l.getSpectrumDifference(10000);
	}
	
	public static void testObjectIO(){
		String filename = "afumigatus_cmp_20.mgf";
		SpectrumLib l = new SpectrumLib(filename, "MGF");
		l.printStat(SpectrumLib.NODETAIL);
		l.removeModSpectra();
		l.removeSingle();
		System.out.println("Original Library");
		l.printStat(SpectrumLib.NODETAIL);
		System.out.println("now writing the library to file");
		l.writeLibToFile("af_lib.obj", l);
		System.out.println("reading back the library from file");
		SpectrumLib l2 = SpectrumLib.readLibFromFile("af_lib.obj");
		l2.printStat(SpectrumLib.NODETAIL);
	}
	
	public static void testCreateMix(){
		String filename = ".\\MSPLib\\Lib\\human.msp";
		//SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		Spectrum mix;
		System.out.println("Original Library: ");
		lib1.removeModSpectra();
		lib1.printStat(SpectrumLib.NODETAIL);
		lib1.toNormVector(1, 0.5, 2000);
		lib1.normIntensity();
		SpectrumLib lib2 = lib1.Divide();
		lib2.printStat(SpectrumLib.NODETAIL);
		//SpectrumLib mixlib = lib2.createMix("mixlibpeptide.dat", 100000, 0.3, 1, 2000); //trying to created the harder set to test various aspect of mix spectrum searching
		//SpectrumLib mixlib = lib2.createMix(10000, 1, 1, 0.3, 1, 2000);
		SpectrumLib mixlib = lib2.createMix("name.txt", 10000, 1, 1, 0.5, 1, 2000);
		System.out.println("The Mixture Library:");
		mixlib.printStat(SpectrumLib.NODETAIL);
		//SpectrumLib.writeLibToFile("mixtureLibrary.obj", mixlib);
		//System.out.println(mixlib); //save it both way just to be sure
		//mixlib.printLib();
	}
	
	public static void testWindFilterPeaks(){
		String filename = ".\\mixture_compressed\\new80min.mgf";
		String fileMix = ".\\mixture_compressed\\new2min.mgf";
        String newfile1 = ".\\mixture_compressed\\out1.mgf";
        String newfile2 = ".\\mixture_compressed\\out2.mgf";
		SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		lib1.printStat(SpectrumLib.NODETAIL);
		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
		mixlib.printStat(SpectrumLib.NODETAIL);
		lib1.windowFilterPeaks(50, 100);
		mixlib.windowFilterPeaks(50, 100);
		lib1.printLibToFile(newfile1, lib1);
		mixlib.printLibToFile(newfile2, mixlib);
	}
	public static void testNumberOfPeaks(){
		String filename = ".\\MSPLib\\Lib\\human.msp";
		//SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		double avg = 0;
		lib1.removeModSpectra();
		lib1.normIntensity();
		//avg = lib1.getPeakCounts(1);
		System.out.println("average is: " + avg);
		System.out.println("removing some peaks");
		lib1.filterPeaks(100);
		lib1.getPeakCounts(1);
	}
	
	public static void testProjectCosine(){
		if(SpectrumLib.DETAIL){
			//return;
		}
		//String filename = "human_cmp_20.mgf";
		//String filename = "afumigatus_cmp_20.mgf";
		String filename = ".\\MSPLib\\Lib\\human.msp";
		//SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		System.out.println("Original Library: ");
		lib1.printStat(SpectrumLib.NODETAIL);
		System.out.println("Removing Spectra with modification");
		lib1.removeModSpectra();
		lib1.printStat(SpectrumLib.NODETAIL);
		System.out.println("==========================================================================================");
		System.out.println();
		System.out.println("Calculating within group similarity:");
		//lib1.getSpectrumSimilarity();
		System.out.println("==========================================================================================");
		System.out.println();
		System.out.println("Calculating btw group similarity/difference:");
		//lib1.getSpectrumDifference(50000);
		System.out.println("==========================================================================================");
		System.out.println();
		System.out.println("Diving Spectrums library into two sets");
		lib1.normIntensity();
		SpectrumLib lib2 = lib1.Divide();
		System.out.println("Dataset 1:");
		lib1.printStat(SpectrumLib.NODETAIL);
		System.out.println("Dataset 2:");
		lib2.printStat(SpectrumLib.NODETAIL);
		SpectrumLib mixlib = lib2.createMix(1000);
		System.out.println("The Mixture Library:");
		mixlib.printStat(SpectrumLib.NODETAIL);
		System.out.println("Calculating projected similarity: ");
		
		Iterator<List<Spectrum>> it = mixlib.getSpectrumLibrary().values().iterator();
		String[] peptides;
		Spectrum mix, candidate;
		List<Spectrum> spects;
		int i = 0;
		while(it.hasNext()){
			mix = it.next().get(0);
			//System.out.println(mix);
			peptides = mix.peptide.split(" & ");
			spects = lib1.getSpectra(peptides[0]);
			for(i = 0; i < spects.size(); i++){
				candidate = spects.get(i).toNormVector(1, 0.5, 2000);
				//System.out.println(candidate);
				//System.out.print("" + peptides[0]);
				System.out.println("" + peptides[0]+ "\t" + candidate.projectedCosine(mix));
			}
			
			spects = lib1.getSpectra(peptides[1]);
			System.out.print("" + peptides[1]);
			for(i = 0; i < spects.size(); i++){
				candidate = spects.get(i).toNormVector(1, 0.5, 2000);
				//System.out.println(candidate);
				//System.out.print("" + peptides[1]);
				System.out.println("" + peptides[1] + "\t" + candidate.projectedCosine(mix));
			}
			//if(peptides[0].equals("DAVTYTEHAK")){
				//return;
			//}
		}
		System.out.println("==========================================================================================");
		System.out.println();
		System.out.println("Calculating Similarity ranks: ");
		int r1, r2;
		lib1.toNormVector(1, 0.5, 2000);
		Iterator mixspects = mixlib.iterator();
		while(mixspects.hasNext()){
			mix = (Spectrum)mixspects.next();
			//System.out.println(mix);
			r1 = lib1.psimilarityRank(mix, mix.peptide);
			System.out.println("" + mix.peptide + ":\t" + r1);		
		}
		
		System.out.println("==========================================================================================");
		System.out.println();
		mixlib = lib2.createMix(50, 1, 1);
		System.out.println("The Mixture Library3:");
		mixlib.printStat(SpectrumLib.NODETAIL);
		System.out.println();
		System.out.println("calculating mix spectrum similarity: ");
		mixspects = mixlib.iterator();
		while(mixspects.hasNext()){
			mix = (Spectrum)mixspects.next();
			//System.out.println(mix);
			peptides = mix.peptide.split(" & ");
			lib1.mixSpectrumSimilarity(mix, peptides[0], peptides[1]);			
		}
	}
	
	public static void testPsimilarity(){
		String filename = ".\\MSPLib\\Lib\\human.msp";
		//String filename = "afumigatus_cmp_20.mgf";
		//SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		Spectrum mix;
		System.out.println("Original Library: ");
		lib1.removeModSpectra();
		lib1.printStat(SpectrumLib.NODETAIL);
		//lib1.normIntensity();
		//lib1.getSpectrumSimilarity(lib1.spectrumList.size());
		//lib1.getSpectrumDifference(50000);
		lib1.toNormVector(1, 0.5, 2000);
		SpectrumLib lib2 = lib1.Divide();
		lib1.printStat(SpectrumLib.NODETAIL);
		lib2.printStat(SpectrumLib.NODETAIL);
		//SpectrumLib mixlib = lib2.createRandomMix(15000, 1, 0.1, 0, 0.0000000001, 10000);
		//SpectrumLib mixlib = lib2.createRandomMix(1000, 1, 0.1, 0.00001, 1, 10000);
		SpectrumLib mixlib = lib2.createMix(".\\mixture_compressed\\testcase10000.name", 1000, 1, 0.1, 0, 1, 2000);
		//SpectrumLib mixlib = lib2.createMix(1000, 1, 1, 0.1, 1, 2000);
		//lib1.normIntensity();
		//lib2.normIntensity();
		//mixlib.normIntensity();
		//mixlib.filterPeaks(250);
		System.out.println("The Mixture Library:");
		mixlib.printStat(SpectrumLib.NODETAIL);
		Iterator mixspects = mixlib.iterator();
		//SpectrumLib mixlib2 = lib2.createMix(50000,1,1,0.2,1);
		//mixlib2.printStat(NODETAIL);
		//Iterator mixspects = mixlib2.iterator();
		System.out.println("Calculating Similarity ranks: ");
		int r1;
		mix = (Spectrum)mixspects.next();
		if(false){
		String[] peps = mix.peptide.split(" & ");
		Spectrum s1 = lib1.getSpectra(peps[0]).get(0);
		Spectrum s2 = lib1.getSpectra(peps[1]).get(0);
		System.out.println("score: " + s1.projectedCosine(mix));
		System.out.println("score: " + s2.projectedCosine(mix));
		}
		if(true){
			while(mixspects.hasNext()){
				mix = (Spectrum)mixspects.next();
				String[] peps = mix.peptide.split(" & ");
				Spectrum s1 = lib1.getSpectra(peps[0]).get(0);
				Spectrum s2 = lib2.getSpectra(peps[0]).get(0);
				Spectrum s3 = lib1.getSpectra(peps[1]).get(0);
				Spectrum s4 = lib2.getSpectra(peps[1]).get(0);
				System.out.println("nonmix similarity: " +  s1.projectedCosine(s2) + "\t" + s3.projectedCosine(s4));
				//lib1.mixSpectrumSimilarity(mix, peps[0], peps[1]);
				//lib1.mixSpectrumDifference(mix, peps[0], peps[1], 1000);
				//System.out.println(mix);
				r1 = lib1.psimilarityRank(mix, mix.peptide);
				//r1 = lib1.psimilarityRank(s2, s2.peptide);
				//r1 = lib1.psimilarityRank(s4, s4.peptide);
				//System.out.println("" + mix.peptide + ":\t" + mix.score + "\t"+ r1);		
			}
		}
	}
	
	public static void testSolution(){
		String filename = ".\\MSPLib\\Lib\\human.msp" ;
		SpectrumLib l = new SpectrumLib(filename, "MSP") ;
		Spectrum mix ;
    		
		l.removeModSpectra() ;
		l.removeSingle() ;
	
		l.toNormVector(1, 0.5, 2000);
		
		SpectrumLib l2 = l.Divide() ;
		//SpectrumLib lib1 = l2.createRandomMix(100, 1, 0.1, 0.00001, 1, 2000) ;
		//l2.filterPeaks(25);
		SpectrumLib lib1 = l2.createMix("bin2.name", 100, 1, 0.1, 0, 1, 2000 ) ;
			
		//l2.normIntensity() ;
		//lib1.normIntensity() ;
		lib1.filterPeaks(125);
		List temp = lib1.getAllSpectrums() ;
		String peptides[];
		Spectrum a, b ;
		//l.filterPeaks(40);
		System.out.println("Start searching");
		for (int j = 0; j < 100; j++) {
			mix = (Spectrum)temp.get(j) ;
			//System.out.println(mix);
			peptides = mix.peptide.split(" & ") ;
			//System.out.print("Solution:\t") ;
			//l.solution(mix) ;
			//System.out.print("Max Solution:\t") ;
			//l.mixSpectrumSimilarity(mix, peptides[0], peptides[1]);
			l.maxSolution(mix);
		}		
	}
	
	public static void testPealingAlg(double ratio, int top, int win){
		String filename = ".\\MSPLib\\Lib\\human.msp" ;
		SpectrumLib l = new SpectrumLib(filename, "MSP") ;
		Spectrum mix ;
    		
		l.removeModSpectra() ;
		l.removeSingle() ;
	
		
		double accuracy = 0;
		SpectrumLib l2 = l.Divide() ;
		//SpectrumLib lib1 = l2.createRandomMix(1000, 1, 0.1, 0.00001, 1, 2000) ;
		SpectrumLib lib1 = l2.createMix("bin2.name", 1000, 1, ratio, 0, 1, 2000) ;			
		int size = lib1.getSpectrumList().size();
		Iterator<Spectrum> itr = l.iterator();
		Spectrum s;
		//l.windowFilterPeaks(top, win);
		//lib1.windowFilterPeaks(top*2, win);
	
		l.toNormVector(1, 0.5, 2000);
		lib1.toNormVector(1, 0.5, 2000);
		
		//l.normIntensity();
		//lib1.normIntensity();
		
		itr = lib1.iterator();
		long start = (new GregorianCalendar()).getTimeInMillis();
		String peptides[];
		Vector<Spectrum> answers;
		Spectrum a, b;
		//lib1.filterPeaks(110);
		//l.filterPeaks(100);
		while(itr.hasNext()){
			mix = itr.next();
			answers = l.peals(mix, 2);
			a = answers.firstElement();
			b = answers.lastElement();
			peptides = mix.peptide.split(" & ");
			if(peptides[0].equals(a.peptide) && peptides[1].equals(b.peptide)){
				accuracy += 1/(double)(size);
				System.out.println(mix.peptide + " answer: " + a.peptide + " & " + b.peptide + " best score: " + a.score +", " + b.score + "\t" + "t");
			}else if(peptides[0].equals(b.peptide) && peptides[1].equals(a.peptide)){
				accuracy += 1/(double)(size);
				System.out.println(mix.peptide + " answer: " + a.peptide + " & " + b.peptide + " best score: " + a.score +", " + b.score + "\t" + "t");
			}else{
				System.out.println(mix.peptide+ " answer: " + a.peptide + " & " + b.peptide + " best score: " + a.score +", " + b.score + "\t" + "f");
			}
		}
		System.out.println("Summary of Result: " + accuracy*100 + "%");
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}

	public static void testSearchAndBoundAlg(){
		//String filename = "afumigatus_cmp_20.mgf";
		String filename = ".\\MSPLib\\Lib\\human.msp";
//		String filename = ".\\mixture_compressed\\new80min.mgf";
//		String fileMix = ".\\mixture_compressed\\new2min.mgf";

//        SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
//		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");

        SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		lib1.removeModSpectra();
		SpectrumLib lib2 = lib1.Divide();
//		SpectrumLib	mixlib = lib2.createRandomMix(1000, 1, 0.1, 0.0001, 1, 1000);
		SpectrumLib	mixlib = lib2.createMix(".\\mixture_compressed\\testcase10000.name", 1000, 1, 1, 0.0001, 1, 1000);
		//SpectrumLib	mixlib = lib2.createMix("bin2.name", 1000, 1, 0.1, 0.0001, 1, 1000);

//		lib1.windowFilterPeaks(15, 50);
//		mixlib.windowFilterPeaks(15, 50);

		lib1.printStat(SpectrumLib.NODETAIL);
		lib1.toNormVector(1, 0.5, 2000);
		mixlib.toNormVector(1, 0.5, 2000);
		lib1.normIntensity();
		mixlib.normIntensity();
		Spectrum temp = null;
		Vector v = null;
		Hashtable <String, Vector <Spectrum>> newTable = new Hashtable();
		int i;
//		for(i = 0; i < 500; i++){
//			temp = lib2.spectrumList.get(i); 
//			lib1.spectrumLibrary.remove(temp.peptide);
//			lib2.spectrumLibrary.remove(temp.peptide);
//			v = new Vector<Spectrum>();
//			v.add(temp);
//			temp.peptide = temp.peptide + " & " + temp.peptide; //conform to mixture pepetide name format
//			newTable.put(temp.peptide, v);
//		}
		//remove half and not remove half
//		for(; i < 500; i++){
//			temp = lib2.spectrumList.get(i); 
//			v = new Vector<Spectrum>();
//			v.add(temp);
//			temp.peptide = temp.peptide + " & " + temp.peptide; //conform to mixture pepetide name format
//			newTable.put(temp.peptide, v);
//		}
//		lib1.spectrumList = lib1.getAllSpectrums();
//		lib2.spectrumList = lib2.getAllSpectrums(); 
//		SpectrumLib mixlib = lib2.createRandomMix(1000, 1, 0.00001, 0.001, 1, 10000);
//		mixlib.spectrumLibrary.putAll(newTable);
//		mixlib.spectrumList = mixlib.getAllSpectrums();
		
//		SpectrumLib mixlib = lib2.createMix("bin2.name", 1000, 1, 0.5, 0.00001, 1, 10000);
//		lib1.windowFilterPeaks(20, 50);
		//System.out.println("mix start: ");
//		mixlib.windowFilterPeaks(20, 50);
		//mixlib.filterPeaks(105);
		//lib1.normIntensity();
		//mixlib.normIntensity();
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
			System.out.println("searching matches to: " + mix.peptide);
			candidate = lib1.searchAndBoundLib2(mix);
//			answer = lib1.findAnswerInSortedSpectrumList(mix, mix.peptide);
			peps = mix.peptide.split(" & ");
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
//			if(answer.score > candidate.score){
//				System.out.println("Do not seems to find the optimal known answer in the DB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
//			}
//			System.out.println(mix);
//			System.out.println("putative	 answer: ");
//			candidate = candidate.toNormVector(1, 0.5, 2000);
//			System.out.println(candidate);
			i++;
			if(i == -1) {
				return;
			}
			
		}
		System.out.println("Summary of Result: " + accuracy*100 + "%");
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
//	public static void test(){
//		String filename = ".\\mixture_compressed\\20070910-A-6mixforMS2_data22.mgf";
//		filename = ".\\mixture_compressed\\new80min.mgf";
//		
//		String fileMix = ".\\mixture_compressed\\20080327-B-Bal6mix18mix_Data07.mgf";
//        //fileMix = ".\\mixture_compressed\\20080328-B-6mixinfuse_Data03.mgf";
//		//String fileMix = filename;
//		SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
//		lib1.printStat(NODETAIL);
//		//lib1.toNormVector(1, 0.5, 2000);
//		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
//		mixlib.printStat(NODETAIL);
//		//mixlib.toNormVector(1, 0.5, 2000);
//		//mixlib = lib1.createRandomMix(100, 1, 1, 0.0000001, 1, 2000);
//		Iterator<Spectrum> it = mixlib.iterator();
//		Spectrum mix, candidate;
//		Vector<Spectrum> candidates;
//		long start = (new GregorianCalendar()).getTimeInMillis();
//		int i = 0;
//		lib1.getConnectedSpectrum(ExperimentalMixtureAnalysis.getSpectrumMixtureGraph(lib1, 2, 0.5, 1000));
//		while(it.hasNext()){
//			mix = it.next();
//			//System.out.println("mix");
//			//System.out.println(mix);
//			//candidate = lib1.searchAndBoundLib(mix);
//			//candidate = lib1.filterAndSearch(mix);
//			//candidates = lib1.peals(mix, 2);
//			//System.out.println(mix.peptide + "\t" + "putative pair: " + candidate.peptide + "\t" + candidate.score + "\t" + candidate.modMass);			
//			//for(i = 0; i < candidates.size(); i++){
//			//	candidate = candidates.get(i);
//			//	System.out.println(candidate.peptide + "\t" + candidate.score + "\t" + candidate.parentMass);
//			//	System.out.println(candidate);
//			//}
//			//System.out.println();
//		}
//		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
//	}


	
	public static void testPealControl(){
		String filename = ".\\MSPLib\\Lib\\human.msp" ;
		SpectrumLib l = new SpectrumLib(filename, "MSP") ;
		Spectrum mix ;
    		
		l.removeModSpectra() ;
		l.removeSingle() ;
	
		l.toNormVector(1, 0.5, 2000);
		double accuracy = 0;
		SpectrumLib l2 = l.Divide() ;
		//SpectrumLib lib1 = l2.createRandomMix(100, 1, 0.1, 0.00001, 1, 2000) ;
		SpectrumLib lib1 = l2.createMix("bin2.name", 10, 1, 1, 0, 1, 2000 ) ;			
		int size = lib1.getSpectrumLibrary().size();
		Iterator<Spectrum> itr = l.iterator();
		Spectrum s;
		
		itr = lib1.iterator();
		String peptides[];
		Vector<Spectrum> answers;
		Spectrum a, b;
		while(itr.hasNext()){
			mix = itr.next();
			answers = l.peals(mix, 2);
			a = answers.firstElement();
			b = answers.lastElement();
			peptides = mix.peptide.split(" & ");
			if(peptides[0].equals(a.peptide) && peptides[1].equals(b.peptide)){
				accuracy += 1/lib1.getSpectrumList().size();
				System.out.println(mix.peptide + " answer: " + a.peptide + " & " + b.peptide + " best score: " + a.score +", " + b.score + "\t" + "t");
			}else if(peptides[0].equals(b.peptide) && peptides[1].equals(a.peptide)){
				accuracy += 1/(size);
				System.out.println(mix.peptide + " answer: " + a.peptide + " & " + b.peptide + " best score: " + a.score +", " + b.score + "\t" + "t");
			}else{
				System.out.println(mix.peptide+ " answer: " + a.peptide + " & " + b.peptide + " best score: " + a.score +", " + b.score + "\t" + "f");
			}
		}
	}
	
	public static void testCosineDev(){
		String filename = ".\\MSPLib\\Lib\\human.msp" ;
		SpectrumLib l = new SpectrumLib(filename, "MSP") ;
		l.removeModSpectra();
		l.toNormVector(1, 0.5, 2000);
		List spectrums = new Vector();
		spectrums.addAll(l.getSpectrumLibrary().values());
		Vector v1, v2;
		double[] stat;
		for(int i = 0; i < spectrums.size(); i++){
			for(int j = 0; j < spectrums.size(); j++){
				v1 = (Vector)spectrums.get(i);
				v2 = (Vector)spectrums.get(j);
				stat = cosineStat(v1, v2);
				System.out.println("mean: " + stat[0] + "\t" + "dev: " + stat[1]);
			}
		}
	}
		
	public static void cosineSpeedTest(){
		//String filename = "afumigatus_cmp_20.mgf";
		String filename = ".\\MSPLib\\Lib\\human.msp";
        SpectrumLib lib1 = new SpectrumLib(filename, "MSP");
		lib1.removeModSpectra();
		lib1.printStat(SpectrumLib.NODETAIL);
		lib1.normIntensity();
		lib1.toNormVector(1, 0.5, 2000);
		long MAXIT = 30000;
		Spectrum s1, s2;
		s1 = lib1.getRandomSpectrum();
		s2 = lib1.getRandomSpectrum();
		System.out.println("size of spectrum: " +  s1.getPeak().size());
		System.out.println("size of spectrum: " + s2.getPeak().size());
		long start = (new GregorianCalendar()).getTimeInMillis();
		double sim;
		for(int i = 0; i < MAXIT; i++){
			sim = s1.cosineSim(s2);
		}
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start) + "msecs");
	}
	
	private static double[] cosineStat(Vector v1, Vector v2){
		double mean = 0, var = 0;
		double[] cosines = new double[v1.size()*v2.size()];
		int k = 0;
		Spectrum s1, s2;
		for(int i = 0; i < v1.size(); i++){
			s1 = (Spectrum)v1.get(i);
			for(int j = 0; j < v2.size(); j++){
				s2 = (Spectrum)v2.get(j);
				cosines[k] = s1.cosineSim(s2);
				mean += cosines[k];
				k++;
			}
		}
		mean = mean / k;
		for(int i = 0; i < cosines.length; i++){
			var += Math.pow(cosines[i]-mean, 2);
		}
		double[] ret = {mean, var};
		return ret;
	}
	
	public static void main(String[] args){
		//testReadMGF();
		// testObjectIO();
		//testProjectCosine();
//		testPsimilarity();
		//testCreateMix();
		testSearchAndBoundAlg();
		//testWeightedCosineSim();
		//testResidualAlpha();
//		testPealingAlg(0.5, 5, 25);
//		testPealingAlg(0.5, 10, 25);
//		testPealingAlg(0.5, 15, 25);
//		testPealingAlg(0.5, 5, 50);
//		testPealingAlg(0.5, 10, 50);
//		testPealingAlg(0.5, 15, 50);
//		testPealingAlg(0.5, 5, 100);
//		testPealingAlg(0.5, 10, 100);
//		testPealingAlg(0.1, 15, 100);
//		testSolution();
		//test();
		//testWindFilterPeaks();
		// SpectrumLib.testReverseSearchMix();
		// testReverseSearchAllPair();
		//SpectrumLib.testRandomRevSearchMix();
		//testReverseSearchControl();
		//testNumberOfPeaks();
		//testSearchNBoundControl();
		//testCosineDev();
		//testExperimentalData();
		//statExperimentalData();
		//cosineSpeedTest();
		
	}
}
