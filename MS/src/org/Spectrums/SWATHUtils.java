package org.Spectrums;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import org.systemsbiology.jrap.stax.Scan;

import IO.MZXMLReader;


/**
 * Some utility function for SWATH data
 * @author Jian Wang
 *
 */
public class SWATHUtils {
	public static final int DALTON = 1;
	public static final int PPM = 2;
	
	
	public static boolean checkPrecursor(double ms2, double precursor, double windowWidth){
		 double diff = precursor - ms2 + 4.0;
		 return diff > 0 && diff < windowWidth;
	}
	
	public static boolean checkMass(double ms2, double precursor, double tolerance, int mode){
		return SpectrumUtil.checkMass(ms2, precursor, tolerance, mode);
	}
	
	
	public static int[] peakIntDistr(Spectrum s, double min, double max, double binWidth){
		List<Peak> sorted = new ArrayList();
		sorted.addAll(s.getPeak());
		Collections.sort(sorted, PeakIntensityComparator.comparator);
		int bins = (int)Math.ceil((max - min)/binWidth);
		int[] counts = new int[bins];
		for(int i = 0; i < sorted.size(); i++){
			
		}
		return null;
	}
	
	public static double getRT(Scan s){
		String rt = s.getHeader().getRetentionTime();
		return Double.parseDouble(rt.substring(2, rt.length()-1));
	}
	
	/**
	 * get the corresponding swath ms1 scan
	 * @param s
	 * @return
	 */
	public static int getSWATHMS1Scan(Spectrum s){
		return (int)(s.scanNumber / 35.0)*35+1;
	}
	
	//convert a spectral library to the peakview ion library text format
	public static void getIonLibrary(String spectrumFile){
		SpectrumLib lib = new SpectrumLib(spectrumFile, "MGF");
		List<Spectrum> specList = lib.getAllSpectrums();
		getIonLibrary(specList);
	}
	
	public static void getIonLibrary(List<Spectrum> specList){
		getIonLibrary(specList, "..//mixture_linked//out.txt");
	}
	
	public static void getIonLibrary(List<Spectrum> specList, String outFile){
		getPeakViewIonLibrary(specList, outFile);
	}
	
	public static void getPeakViewIonLibrary(List<Spectrum> specList, String outFile){
		int counter = 0;
		int ProteinCounter = 1; //not sure if this is protein group number but in ionlibrary seems unique to each proteins
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		String header = "Q1\tQ3\tRT_detected\tprotein_name\tisotype\trelative_intensity\tstripped_sequence\tmodification_sequence\tprec_z\tfrg_type\tfrg_z\tfrg_nr\tiRT\tuniprot_id\tscore\tdecoy\tprec_y\tconfidence\tshared\tN\tmods\tnterm\tcterm";
		out.println(header);
		Map<String, Integer> proteinMap = new HashMap<String, Integer>();		
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			//if(s.peptide.contains("+") || s.spectrumName.contains("DECOY")){
			//if(s.spectrumName.contains("DECOY") || 
			//		s.peptide.matches(".*[^C]\\[.*")){
			if(s.peptide.contains("+")){
				//System.out.println("Matching: " + s.peptide + "\t" + s.spectrumName);
				continue;
			}
			SimpleMatchingGraph g = SpectrumUtil.getMatchGraph(s, 0.05);
			TreeSet sortedSet = new TreeSet(PeakIntensityComparator.comparator);
			sortedSet.addAll(g.vertexSet(SimpleMatchingGraph.Observed));
			int matchCount = 0;
			StringBuffer outbuffer = new StringBuffer();
			//for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			//s.rt = (s.rt+60)*72;
			String stripped = s.peptide.replaceAll("[0-9\\.\\+\\-\\[\\]]", "");
			String modSeq = stripped.replaceAll("C", "C[CAM]");
			//s.protein = s.protein.replaceAll("/", "");
			for(Iterator<LabelledPeak> it = sortedSet.descendingIterator(); it.hasNext();){
				Peak p = it.next();
				List neighs = g.getNeighbors(p);
				LabelledPeak closestP = null;
				double smallestErr = 10000;
				for(int n = 0; n < neighs.size(); n++){
					LabelledPeak currPeak = (LabelledPeak)neighs.get(n);
					double currDiff = Math.abs(currPeak.getMass()-p.getMass());
					closestP = currDiff < smallestErr ? currPeak : closestP;
					smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
				}
				int protInd = 1;
				
				String[] tokens = s.spectrumName.split("\\s+");
//				if(s.rt == Double.MIN_VALUE){
//					s.rt = Double.parseDouble(tokens[5].substring(2, tokens[5].length()-1));
//				}
				if(s.protein == null){
					s.protein = tokens[7];	
				}	
				if(proteinMap.containsKey(s.protein)){
					protInd = proteinMap.get(s.protein);
				}else{
					protInd = ProteinCounter++;
					proteinMap.put(s.protein, protInd);
				}
				double Int = 10000;
				if(closestP != null 
						//&& (closestP.getType().equals("b") || closestP.getType().equals("y"))
						&&		closestP.getPos() > 1){
					matchCount++;
					outbuffer.append(s.parentMass + "\t");
					outbuffer.append(p.getMass() + "\t");
					outbuffer.append(s.rt/60 + "\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(" " + "\t");
					outbuffer.append(p.getIntensity() + "\t");
					outbuffer.append(stripped +"\t");					
					outbuffer.append(modSeq +"\t");
					//System.out.print(s.peptide + "\t");
					outbuffer.append(s.charge  + "\t");
					outbuffer.append(closestP.getType() + "\t");
					outbuffer.append(closestP.getCharge() + "\t");
					outbuffer.append(closestP.getPos() + "\t");
					outbuffer.append(s.rt/60 + "\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(-1 + "\t");
					outbuffer.append("False" + "\t");
					outbuffer.append(Int + "\t");
					outbuffer.append("1" + "\t");
//					if(s.score == 1)
//						outbuffer.append("False" + "\t");
//					else
//						outbuffer.append("True" + "\t");
					outbuffer.append("False\t");
					outbuffer.append(protInd + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("" + "\t");
					outbuffer.append("\n");
				}
				
				if(matchCount > 500){
					break;
				}

			}
			if(matchCount >= 5){
				out.print(outbuffer);
			}else{
				//System.out.println("warning less than 10 peaks in library spectrum");
			}
			counter++;
		}
		out.flush();
		out.close();
	}
	
	
	public static void getOpenSwathAssayLibrary(List<Spectrum> specList, String outFile){
		int counter = 0;
		int ProteinCounter = 1; //not sure if this is protein group number but in ionlibrary seems unique to each proteins
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		String header = "PrecursorMz\tProductMz\tTr_recalibrated\ttransition_name\tCE\tLibraryIntensity\ttransition_group_id\tdecoy\tPeptideSequence\tProteinName\tAnnotation\tFullUniModPeptideName\tMissedCleavages\tReplicates\tNrModifications\tPrecursorCharge\tGroupLabel";
		out.println(header);
		Map<String, Integer> proteinMap = new HashMap<String, Integer>();		
		int Ind = 0;
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			//System.out.println("Number of peaks: " + s.getPeak().size());
			s.mergePeaks(s, 0.1);
			s.removePeaksInMass(310, 2500);
			//s.rt = s.rt/60;
			//System.out.println("RT " + s.rt);
			if(s.peptide.contains("+")){
			//if(s.spectrumName.contains("DECOY")){
				continue;
			}
			SimpleMatchingGraph g = SpectrumUtil.getMatchGraph(s, 0.05);
			TreeSet sortedSet = new TreeSet(PeakIntensityComparator.comparator);
			sortedSet.addAll(g.vertexSet(SimpleMatchingGraph.Observed));
			int matchCount = 0;
			StringBuffer outbuffer = new StringBuffer();
			//for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Ind++;
			for(Iterator<LabelledPeak> it = sortedSet.descendingIterator(); it.hasNext();){
				Peak p = it.next();
				List neighs = g.getNeighbors(p);
				LabelledPeak closestP = null;
				double smallestErr = 10000;
				for(int n = 0; n < neighs.size(); n++){
					LabelledPeak currPeak = (LabelledPeak)neighs.get(n);
					double currDiff = Math.abs(currPeak.getMass()-p.getMass());
					closestP = currDiff < smallestErr ? currPeak : closestP;
					smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
				}
				int protInd = 1;
				String[] tokens = s.spectrumName.split("\\s+");
			
				
				if(s.rt == Double.MIN_VALUE){
					s.rt = Double.parseDouble(tokens[5].substring(2, tokens[5].length()-1));
					s.rt = s.rt/60;
				}
				
				if(s.protein == null){
					s.protein = tokens[7];	
				}
				s.protein = tokens[7];
				if(proteinMap.containsKey(s.protein)){
					protInd = proteinMap.get(s.protein);
				}else{
					protInd = ProteinCounter++;
					proteinMap.put(s.protein, protInd);
				}
				double Int = 10000;
				String pep = Utils.StringUtils.getStrippedSeq(s.peptide);
				Peptide peptide = new Peptide(s.peptide, s.charge);
				if(closestP != null 
						&& (closestP.getType().equals("b") || closestP.getType().equals("y"))
						&&		closestP.getPos() > 1){
					matchCount++;
					//outbuffer.append(s.parentMass + "\t");
					outbuffer.append(peptide.getParentmass() + "\t");
					outbuffer.append(p.getMass() + "\t");
					outbuffer.append(s.rt + "\t");
					outbuffer.append(Ind+"_"+s.peptide+"/"+s.charge+"_"+closestP.getType()+closestP.getPos()+"/"+closestP.getCharge() +"\t" );
					outbuffer.append("27"+"\t");
					outbuffer.append(p.getIntensity()+"\t");
					outbuffer.append(Ind+"_"+s.peptide+"/"+s.charge+"\t");
					if(s.spectrumName.contains("DECOY")){
						outbuffer.append("1"+"\t");
					}else{
						outbuffer.append("0"+"\t");
					}
					outbuffer.append(pep +"\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(closestP.getType()+closestP.getPos() + "\t");
					if(s.peptide.contains("C")){
						StringBuffer modPep = new StringBuffer(pep);
						int index = modPep.indexOf("C");
						//System.out.println("ModPep " + modPep);
						//System.out.println("Index1: " + index);						
						while(index > -1){
							modPep.insert(index+1, "(UniMod:4)");
							//System.out.println("Index: " + index);
							//System.out.println("ModPep:" + modPep);
							index = modPep.indexOf("C", index+1);
//							System.out.println("Index: " + index);
							
						}
						outbuffer.append(modPep.toString() + "\t");
					}else{
						outbuffer.append(pep + "\t");
					}
					outbuffer.append("0" + "\t");
					outbuffer.append("0" + "\t");
					outbuffer.append("0" + "\t");
					outbuffer.append(s.charge + "\t");
					outbuffer.append("light" + "\n");
				}
				if(matchCount > 5){
					break;
				}

			}
			if(matchCount >= 5){
				out.print(outbuffer);
			}else{
				//System.out.println("warning less than 10 peaks in library spectrum");
			}
			counter++;
		}
		out.flush();
		out.close();
	}
	
	//text-based assay library format for skyline
	public static void getSkylineAssayLibrary(List<Spectrum> specList, String outFile){
		int counter = 0;
		int ProteinCounter = 1; //not sure if this is protein group number but in ionlibrary seems unique to each proteins
		PrintStream out = Utils.FileIOUtils.getOutStream(outFile);
		String header = "PrecursorMz\tProductMz\tTr_recalibrated\ttransition_name\tCE\tLibraryIntensity\ttransition_group_id\tdecoy\tPeptideSequence\tProteinName\tAnnotation\tFullUniModPeptideName";
		out.println(header);
		Map<String, Integer> proteinMap = new HashMap<String, Integer>();		
		int Ind = 0;
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			//System.out.println("Number of peaks: " + s.getPeak().size());
			s.mergePeaks(s, 0.1);
			s.removePeaksInMass(310, 2500);
			//s.rt = s.rt/60;
			//System.out.println("RT " + s.rt);
			if(s.peptide.contains("+")){
			//if(s.spectrumName.contains("DECOY")){
				continue;
			}
			SimpleMatchingGraph g = SpectrumUtil.getMatchGraph(s, 0.05);
			TreeSet sortedSet = new TreeSet(PeakIntensityComparator.comparator);
			sortedSet.addAll(g.vertexSet(SimpleMatchingGraph.Observed));
			int matchCount = 0;
			StringBuffer outbuffer = new StringBuffer();
			//for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Ind++;
			for(Iterator<LabelledPeak> it = sortedSet.descendingIterator(); it.hasNext();){
				Peak p = it.next();
				List neighs = g.getNeighbors(p);
				LabelledPeak closestP = null;
				double smallestErr = 10000;
				for(int n = 0; n < neighs.size(); n++){
					LabelledPeak currPeak = (LabelledPeak)neighs.get(n);
					double currDiff = Math.abs(currPeak.getMass()-p.getMass());
					closestP = currDiff < smallestErr ? currPeak : closestP;
					smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
				}
				int protInd = 1;
				String[] tokens = s.spectrumName.split("\\s+");
			
				
				if(s.rt == Double.MIN_VALUE){
					s.rt = Double.parseDouble(tokens[5].substring(2, tokens[5].length()-1));
					//s.rt = s.rt/60;
				}
				
				if(s.protein == null){
					s.protein = tokens[7];	
				}
				s.protein = tokens[7];
				if(proteinMap.containsKey(s.protein)){
					protInd = proteinMap.get(s.protein);
				}else{
					protInd = ProteinCounter++;
					proteinMap.put(s.protein, protInd);
				}
				double Int = 10000;
				String pep = Utils.StringUtils.getStrippedSeq(s.peptide);
				Peptide peptide = new Peptide(s.peptide, s.charge);
				if(closestP != null 
						&& (closestP.getType().equals("b") || closestP.getType().equals("y"))
						&&		closestP.getPos() > 1){
					matchCount++;
					//outbuffer.append(s.parentMass + "\t");
					outbuffer.append(peptide.getParentmass() + "\t");
					outbuffer.append(p.getMass() + "\t");
					outbuffer.append(s.rt + "\t");
					outbuffer.append(Ind+"_"+s.peptide+"/"+s.charge+"_"+closestP.getType()+closestP.getPos()+"/"+closestP.getCharge() +"\t" );
					outbuffer.append("27"+"\t");
					outbuffer.append(p.getIntensity()+"\t");
					outbuffer.append(Ind+"_"+s.peptide+"/"+s.charge+"\t");
					if(s.spectrumName.contains("DECOY")){
						outbuffer.append("1"+"\t");
					}else{
						outbuffer.append("0"+"\t");
					}
					outbuffer.append(pep +"\t");
					outbuffer.append(s.protein + "\t");
					outbuffer.append(closestP.getType()+closestP.getPos() + "\t");
					if(s.peptide.contains("C")){
						StringBuffer modPep = new StringBuffer(pep);
						int index = modPep.indexOf("C");
						//System.out.println("ModPep " + modPep);
						//System.out.println("Index1: " + index);						
						while(index > -1){
							modPep.insert(index+1, "(UniMod:4)");
							//System.out.println("Index: " + index);
							//System.out.println("ModPep:" + modPep);
							index = modPep.indexOf("C", index+1);
//							System.out.println("Index: " + index);
							
						}
						outbuffer.append(modPep.toString() + "\t");
					}else{
						outbuffer.append(pep + "\t");
					}
				}
				if(matchCount > 5){
					break;
				}

			}
			if(matchCount >= 5){
				out.print(outbuffer);
			}else{
				//System.out.println("warning less than 10 peaks in library spectrum");
			}
			counter++;
		}
		out.flush();
		out.close();
	}
	
	
	/**
	 * Generate a ion library for targeted extraction based on identification from multiplexed 
	 * search tools such as M-SPLIT, the identified peptides are then used to do targeted extraction from swath
	 * @param libFile
	 * @param SpectrumFile
	 * @param resultsFile
	 */
	
	public static void getIonLibrary(String libFile, String spectrumFile, String resultsFile, String fastaFile, String outFile){
		
	}
	
	public static void getIonLibrary(String libFile, String spectrumFile, String resultsFile, String fastaFile, String outFile, String outFormat){
		SpectrumLib lib = new SpectrumLib(libFile, "MGF");
		ProteinIDExtractor protID = new ProteinIDExtractor(lib.getAllSpectrums(), fastaFile);
		Map<String, List<String>> pMap = protID.peptideMap;
		//MZXMLReader reader = new MZXMLReader(spectrumFile);
		//TDAStat searchResults = new TDAStat(resultsFile, 1, 4, 6, 8, 32, -1);
		TDAStat searchResults = new TDAStat(resultsFile, 1, 16, 4, 6, 8, 32, -1, false, 0.01, fastaFile);
		List<String> results = Utils.FileIOUtils.createListFromFile(resultsFile);
		List<Spectrum> specList = new ArrayList<Spectrum>();
		int counter = 0;
		Collection<AnnotatedSpectrum> peps = searchResults.getPeptideResults();
		for(Iterator<AnnotatedSpectrum> it = peps.iterator(); it.hasNext();){
			AnnotatedSpectrum s = it.next();
			String pepID = s.peptide + "." + s.charge;
			//System.out.println("peptide: " + pepID);
			if(lib.getSpectra(pepID) != null){
				Spectrum libSpect = lib.getSpectra(pepID).get(0);
				//Spectrum swathSpect = reader.getSpectrum(s.scanNumber);
				//MixSSM mSSM = searchResult.getMatches(swathSpect.scanNumber);
				//if(mSSM == null){
				//	continue;
				//}
				//Spectrum unique = mSSM.getUniqueSpect(libSpect);
				//System.out.println("Spectrum size: " + libSpect.getPeak().size() + "\t" + unique.getPeak().size());
				libSpect.rt = s.rt;// swathSpect.rt;
				libSpect.scaleSpectrum(1/libSpect.sumMagnitude());
				//libSpect.setPeaks(unique.getPeak());
				//String[] tokens = libSpect.spectrumName.split("\\s+");
				//libSpect.rt = Double.parseDouble(tokens[5].substring(2, tokens[5].length()-1));		
				if(pMap.containsKey(libSpect.peptide)){
					List<String> prots = pMap.get(libSpect.peptide);
					libSpect.score = prots.size();
					libSpect.protein = prots.get(0);
					specList.add(libSpect);
					counter++;
				}
				if(counter % 1000 == 0){
					System.out.println("got spec: " + counter);
				}
			}
		}
		System.out.println("Getting library for peptides: " + specList.size());
		if(outFormat.equals("PeakView")){
			getPeakViewIonLibrary(specList, outFile);
			//getIonLibrary(specList, outFile);
		}else if(outFormat.equals("OpenSwath")){
			getOpenSwathAssayLibrary(specList, outFile);
		}else if(outFormat.equals("Skyline")){
				getSkylineAssayLibrary(specList, outFile);
		}else{
			System.err.println("Unknown assay library format");
		}
	}
	

	
	/**
	 * Generate ion library base on fragments ions from swath directly
	 * (as oppose to those from reference library)
	 * @param libFile
	 * @param spectrumFile
	 * @param resultsFile
	 */
	public static void getIonLibraryFromSWATH(String libFile, String spectrumFile, String resultsFile){
		SpectrumLib lib = new SpectrumLib(libFile, "MGF");
		ProteinIDExtractor protID = new ProteinIDExtractor(lib.getAllSpectrums(), "../mixture_linked/database/UPS_plusHuman_plusDecoy.fasta");
		Map<String, List<String>> pMap = protID.peptideMap;
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		List<String> results = Utils.FileIOUtils.createListFromFile(resultsFile);
		List<Spectrum> specList = new ArrayList<Spectrum>();
		int counter = 0;
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String[] result = it.next().split("\\t");
			String ID = result[4]+"."+result[6];
			if(lib.getSpectra(ID) != null){
				Spectrum libSpect = lib.getSpectra(ID).get(0);
				Spectrum swathSpect = reader.getSpectrum(Integer.parseInt(result[1]));
				//swathSpect.windowFilterPeaks2(15, 25);
				//swathSpect.mergePeaks(swathSpect, 0.05);
				Spectrum projection = swathSpect.project(libSpect, 0.05);
				projection.mergePeaks(projection, 0.05);
				//projection = projection.project(libSpect, 0.05);
				projection.rt = swathSpect.rt;
				projection.parentMass = libSpect.parentMass;
				projection.charge = libSpect.charge;
				projection.peptide = libSpect.peptide;
				//System.out.println(libSpect.peptide);
				if(pMap.containsKey(libSpect.peptide)){
					List<String> prots = pMap.get(libSpect.peptide);
					projection.score = prots.size();
					projection.protein = prots.get(0);
					specList.add(projection);
					counter++;
				}
				if(counter % 1000 == 0){
					System.out.println("got spec: " + counter);
				}
			}
		}
		//System.out.println(specList.size());
		getIonLibrary(specList);
	}
	
	
	public static void testGenerateIonLibrary(String specLibFile, String resultFile, String fastaFile, String outFile){
		specLibFile = "../mixture_linked/ACG_swathdevelopment_UPSEcoli_REP234_IDA.mgf";
		resultFile = "../mixture_linked/t0";
		fastaFile = "../mixture_linked/database/UPS_plusEcoli_plusDecoy.fasta";
		outFile = "..//mixture_linked///test_openSWATH/testSkylineLib_withcorrectedPM.txt";
		SpectrumLib lib = new SpectrumLib(specLibFile, "MGF");
		//String ssmResultFile = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPSEcoli_libREp234_msplit_7_pep.txt";
		//getIonLibrary(specLibFile);
		//getIonLibrary(specLibFile, "", resultFile, fastaFile, outFile);
		getOpenSwathAssayLibrary(lib.getAllSpectrums(), outFile);
		//getIonLibraryFromSWATH(specLibFile, spectrumFile, resultFile);
	}
	
	public static void getFilteredSWATHFile(String swathFile, String outFile){
		//swathFile = "../../Downloads/ACG_Nuno-selected/New_Swath_3mu/SWATH_3amu_01.mzXML";
		//outFile = "../mixture_linked/out.mgf";
		MZXMLReader reader = new MZXMLReader(swathFile);
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			while(reader.hasNext()){
				Spectrum s = reader.next();
				s.windowFilterPeaks(15, 25);
				//s.computePeaksZScore(0.5);
				//s.filterPeaksByRankScore(200);
				if(s.getPeak().size() > 10){
					//not printing ranks
					for(int i = 0; i < s.getPeak().size(); i++){
						s.getPeak().get(i).setRank(0);
					}
					out.write(s.toString());
					out.write("\n");
				}
			}
			out.flush();
			out.close();
		}catch(IOException ioe){
			
		}
	}
	
	public static void computeSimLibrary(){
		String libTarget = "../mixture_linked/APSWATH_AP3REPlib_IDA_1pepFDR_plusDecoy2.mgf";
		String libDecoy = "../mixture_linked/ACG_swathdevelopment_UPSEcoli_REP234_IDA_plusDecoy2.mgf";
		SpectrumLib lib = new SpectrumLib(libTarget, "MGF");
		SpectrumLib lib2 = new SpectrumLib(libDecoy, "MGF");
		List<Spectrum> specList = lib.getAllSpectrums();
		List<Spectrum> specList2 = lib2.getAllSpectrums();
		for(int i = 0; i < specList.size(); i++){
			Spectrum s1 = specList.get(i);
			s1.mergePeaks(s1, 0.05);
			s1.sqrtSpectrum();
			for(int j = 0; j < specList.size(); j++){
				Spectrum s2 = specList2.get(j);
				s2.mergePeaks(s2, 0.05);
				s2.sqrtSpectrum();
				if(Math.abs(s1.parentMass - s2.parentMass) < 30){
					double score12 = s1.projectedCosine(s2, 0.05);
					double score21 = s2.projectedCosine(s1, 0.05);
					if(score12 > 0.5 || score21 > 0.5){
						System.out.println(s1.spectrumName + "\t"   + s1.peptide + "\t" + s1.charge + "\tand\t" 
								+ s2.spectrumName + "\t" + s2.peptide + "\t"  + s2.charge + "\tsim:\t" 
								+ score12 + "\t" + score21);
					}
				}
			}
		}
		
	}
	
	public static void main(String[] args){
		testGenerateIonLibrary(args[1], args[2], args[3], args[4]);
		//computeSimLibrary();
		//getFilteredSWATHFile(args[0], args[1]);
		//System.out.println("Unknown command: " + args[0]);
	}
}

