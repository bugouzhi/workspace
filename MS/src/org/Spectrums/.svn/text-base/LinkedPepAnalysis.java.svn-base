package org.Spectrums;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.ArrayList;

/**
 * Containers of some methods that analyze linked peptide
 * @author jian wang
 *
 */
public class LinkedPepAnalysis {
		private static double parentMassTolerance = 0.01;
		public static Vector findLinkCandidateByMass(Spectrum linked, SpectrumLib unlinkedList, double linkerMass, double dangleMass){
			Vector results = new Vector();
			Spectrum s1, s2;
			List<Spectrum> list = unlinkedList.getAllSpectrums();
			for(int i = 0; i < list.size(); i++){
					s1 = list.get(i);
				for(int j = i+1; j < list.size(); j++){
					s2 = list.get(j);
					if(isLinkedCandidate(linked, s1, s2, linkerMass, dangleMass)){
						Spectrum[] pair = {s1, s2};
						results.add(pair);
					}
				}	
			}
			return results;
		}
		
		public static Vector findLinkCandidateByMass(Spectrum linked, List<Peptide> pepList, double linkerMass, double dangleMass){
			Vector results = new Vector();
			Peptide p1, p2;
			for(int i = 0; i < pepList.size(); i++){
					p1 = pepList.get(i);
				for(int j = i+1; j < pepList.size(); j++){
					p2 = pepList.get(j);
					if(isLinkedCandidate(linked, p1, p2, linkerMass, dangleMass)){
						Peptide[] pair = {p1, p2};
						results.add(pair);
					}
				}	
			}
			return results;
		}
		
		private static boolean isLinkedCandidate(Spectrum linked, Spectrum unlinked1, Spectrum unlinked2, double linkerMass, double dangleMass){
			//System.out.println("linked pm:" + linked.parentMass);
			//System.out.println("unliked pm: " + unlinked1.parentMass + "\t" + unlinked2.parentMass);
			boolean massMatch = Math.abs((linked.parentMass*linked.charge - linkerMass)  
						- (unlinked1.parentMass*unlinked1.charge + unlinked2.parentMass*unlinked2.charge - 2*dangleMass)) < parentMassTolerance;
//			if(massMatch){
//				System.out.println("linked charge: " + linked.charge);
//				System.out.println("unlinked charge: " + unlinked1.charge + "\t" + unlinked2.charge);
//			}
			boolean chargeMatch = true && linked.charge == unlinked1.charge + unlinked2.charge;
			boolean annotated = !unlinked1.peptide.contains("Scan") && !unlinked2.peptide.contains("Scan");
			annotated = annotated && unlinked1.peptide.contains("K") && unlinked2.peptide.contains("K");
			return (massMatch & chargeMatch & annotated);
		}
		
		private static boolean isLinkedCandidate(Spectrum linked, Peptide unlinked1, Peptide unlinked2, double linkerMass, double dangleMass){
			//System.out.println("linked pm:" + linked.parentMass);
			//System.out.println("unliked pm: " + unlinked1.parentMass + "\t" + unlinked2.parentMass);
			double target = linked.parentMass*linked.charge;
			double candidate = unlinked1.getParentmass()*unlinked1.getCharge() + unlinked2.getParentmass()*unlinked2.getCharge() 
						+ linkerMass + (linked.charge - unlinked1.getCharge() - unlinked2.getCharge())*Mass.PROTON_MASS
						- 2*dangleMass;
			double massDiff = Math.abs(target - candidate);
			boolean massMatch = massDiff < parentMassTolerance;
//				Math.abs((linked.parentMass*linked.charge - linkerMass)  
//						- (unlinked1.getParentmass()*unlinked1.getCharge() + unlinked2.getParentmass()*unlinked2.getCharge()
//								+ (linked.charge - unlinked1.getCharge() - unlinked2.getCharge())*Mass.PROTON_MASS
//								- 2*dangleMass)) < parentMassTolerance;
			return massMatch;
		}
		
		public static void candidateLinkedMatch(String linked, String unlinked){
			String  filename = unlinked;
			String fileMix = linked;
			String annotationFile =".\\mixture_linked\\trps\\result.txt";
			SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
			lib1.printStat(SpectrumLib.NODETAIL);
			lib1.annoateSpectrumFromInspectFile(annotationFile);
			//lib1.toNormVector(1, 0.5, 2000);
			SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
			mixlib.printStat(SpectrumLib.NODETAIL);
			
			Iterator<Spectrum> it = mixlib.iterator();
			Spectrum mix;
			Spectrum[] candidate;
			Vector candidates;
			long start = (new GregorianCalendar()).getTimeInMillis();
			while(it.hasNext()){
				mix = it.next();
				candidates = findLinkCandidateByMass(mix, lib1, 173.9809, 145.0197);
				for(int i = 0; i < candidates.size(); i++){
					candidate = (Spectrum[])candidates.get(i);
					System.out.println(mix.spectrumName + "\t" + mix.parentMass + "\t" +  mix.charge
							+ "\t" + candidate[0].spectrumName + "\t" + candidate[0].parentMass + "\t" + candidate[0].charge 
							+ "\t" + candidate[1].spectrumName + "\t" + candidate[1].parentMass + "\t" + candidate[1].charge
							+ "\t" + candidate[0].peptide + "\t" + candidate[1].peptide);
				}
			}
			System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		}
		
		public static void candidateLinkedPeptides(String linked, String peptideFile){
			String fileMix = linked;
			List<Peptide> lib1 = loadPeptides(peptideFile);
			//lib1.toNormVector(1, 0.5, 2000);
			SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
			mixlib.printStat(SpectrumLib.NODETAIL);	
			Iterator<Spectrum> it = mixlib.iterator();
			Spectrum mix;
			Peptide[] candidate;
			Vector candidates;
			long start = (new GregorianCalendar()).getTimeInMillis();
			while(it.hasNext()){
				mix = it.next();
				if(mix.charge >= 3){
					candidates = findLinkCandidateByMass(mix, lib1, 173.9809, 0);
					for(int i = 0; i < candidates.size(); i++){
						candidate = (Peptide[])candidates.get(i);
						System.out.println(mix.spectrumName + "\t" + mix.parentMass + "\t" +  mix.charge
							+ "\t" + candidate[0].getPeptide() + "\t" + candidate[0].getParentmass() 
							+ "\t" + candidate[0].getCharge() 
							+ "\t" + candidate[1].getPeptide() + "\t" + candidate[1].getParentmass() 
							+ "\t" + candidate[1].getCharge());
					}
				}
			}
			System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		}
		
		public static List<Peptide> loadPeptides(String peptides){
			List<Peptide> pList = null;
			try{
				BufferedReader bf = new BufferedReader(new FileReader(peptides));
				String currentLine = bf.readLine();
				StringBuffer protein = null;
				String header = null;
				pList = new ArrayList();
				while(currentLine != null){
					Peptide p = new Peptide(currentLine+".1");
					pList.add(p);
					currentLine = bf.readLine();
				}
				bf.close();
				
			}catch(IOException ioe){
				System.out.println(ioe.getMessage());
				System.out.println(ioe.getCause());
			}
			return pList;
		}
		
//		public static void createTargetedSubSet(){
//			//String unlinked = ".\\mixture_linked\\tk090204_WKarzai_DTT_Chymo.mgf";
//			//String unlinked = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
//			//String linked = ".\\mixture_linked\\tk090204_WKarzai_Tryp.mgf";
//			String linked = ".\\mixture_linked\\tk090204_WKarzai_Chymo.mgf";
//			String names = ".\\mixture_linked\\name.txt";
//			String outfile =".\\mixture_linked\\out.txt";
////			SpectrumLib lib1 = new SpectrumLib(unlinked, "MGF");
//			SpectrumLib lib1 = new SpectrumLib(linked, "MGF");
//			SpectrumLib smallLib = lib1.getSubLib(names);
//			smallLib.printLibToFile(outfile, smallLib);
//		}
		
		public static void filterNonLinkedSpectrum(){
//			String spectrumFile = ".\\mixture_linked\\tk090204_WKarzai_DTT_Tryp.mgf";
//			String mixtureFile = ".\\mixture_linked\\tk090204_WKarzai_Tryp.mgf";
			String annotationFile =".\\mixture_linked\\trps\\result.txt";
			String tripletFile =".\\mixture_linked\\triplet.txt";
			
			String spectrumFile = ".\\mixture_linked\\tk090204_WKarzai_DTT_Chymo.mgf";
			String mixtureFile = ".\\mixture_linked\\tk090204_WKarzai_Chymo.mgf";
			
			SpectrumLib mix = new SpectrumLib(spectrumFile, "MGF");
			System.out.println("Done loading spectrum library");
			//lib1.annoateSpectrumFromInspectFile(annotationFile);
			SpectrumLib lib1 = new SpectrumLib(mixtureFile, "MGF");
			System.out.println("Done loading linked-mixture library");
			lib1.toNormVector(1, 0.5, 2000);
			mix.toNormVector(1, 0.5, 2000);
			lib1.normIntensity();
			mix.normIntensity();
			
			Iterator<Spectrum> it = mix.iterator();
			Spectrum m, candidate, answer;
			Vector<Spectrum> spects;
			String[] peps, putativepeps;
			long start = (new GregorianCalendar()).getTimeInMillis();
			double accuracy = 0, size = mix.getAllSpectrums().size();
			Vector<Spectrum> candidates;
			while(it.hasNext()){
				m = it.next();
				if(m.getPeak().size() == 0){ //make sure the mixture has some peaks otherwise it will take long time to search
					continue;
				}
				System.out.println("searching matches to: " + m.peptide);
				candidate = lib1.topSpectrum(m);
				System.out.println("best answer: " + m.peptide + "\t" + m.parentMass + "\t" + m.charge + "\t" + candidate.score
						+ "\t" + candidate.peptide + "\t" + candidate.parentMass + "\t" + candidate.charge);
			}
		}
		
		public static void spectrumLibStat(){
			//String file = ".\\mixture_linked\\Tryp_filtered_linked.mgf";
			String file = ".\\mixture_linked\\tk090204_WKarzai_Tryp.mgf";
			SpectrumLib mix = new SpectrumLib(file, "MGF");
			Iterator<Spectrum> it = mix.iterator();
			Spectrum m;
			while(it.hasNext()){
				m = it.next();
				System.out.println(m.spectrumName + " parentmass: " + m.parentMass + " charge: " + m.charge);
			}
		}
		
		
		public static void main(String[] args){
//			createTargetedSubSet();
//			candidateLinkedMatch(".\\mixture_linked\\tk090204_WKarzai_Chymo.mgf", 
//					".\\mixture_linked\\tk090204_WKarzai_DTT_Chymo.mgf");
//			candidateLinkedMatch(".\\mixture_linked\\tk090204_WKarzai_Tryp.mgf", 
//			".\\mixture_linked\\Try_filtered.mgf");
//			candidateLinkedMatch(".\\mixture_linked\\Tryp_filtered_linked.mgf", 
//			".\\mixture_linked\\tk090204_WKarzai_DTT_Chymo.mgf");
//			filterNonLinkedSpectrum();
			String peptideFile ="..\\mixture_linked\\TrypDigested_filteredPeptides.txt";
			//String peptideFile ="..\\mixture_linked\\Ecoli_allpeptides.txt";
			candidateLinkedPeptides("..\\mixture_linked\\Tryp_filtered_linked.mgf", peptideFile);
//			spectrumLibStat();
		}
}
