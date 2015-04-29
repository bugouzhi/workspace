package org.Spectrums;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import IO.MZXMLReader;

/**
 * This class contain method to analyze multiple SWATH/DIA runs
 * @author Jian
 *
 */
public class MultiSWATHAnalyzer {
	List<String> swathFiles;
	List<MZXMLReader> readers;
	List<TDAStat> results;
	String libFile;
	SpectrumLib lib;
	
	
	public MultiSWATHAnalyzer (String swath1, String swath2, String result1, String result2, String libFile){
		this.swathFiles = new ArrayList<String>();
		this.swathFiles.add(swath1);
		this.swathFiles.add(swath2);
		this.readers = new ArrayList<MZXMLReader>();
		this.readers.add(new MZXMLReader(swath1));
		this.readers.add(new MZXMLReader(swath2));
		this.results = new ArrayList<TDAStat>();
		this.results.add(new TDAStat(result1, 1, 4, 6, 8,31,-1, true));
		this.results.add(new TDAStat(result2, 1, 4, 6, 8,31,-1, true));
		this.libFile = libFile;
		this.lib = new SpectrumLib(libFile, "MGF");
	}
	
	public MultiSWATHAnalyzer(String[] resultsFiles){
		this.results = new ArrayList<TDAStat>(resultsFiles.length);
		int total = resultsFiles.length;
		for(int i = 0; i < total; i++){
			results.add(new TDAStat(resultsFiles[i], 1,4,6,8,31,-1,true));
		}
		
		
	}
	
	public MultiSWATHAnalyzer(String[] resultsFiles, int[][] index, String fasta){
		this.results = new ArrayList<TDAStat>(resultsFiles.length);
		int total = resultsFiles.length;
		boolean useCharge = false;
		for(int i = 0; i < total; i++){
			results.add(new TDAStat(resultsFiles[i], index[i][0], index[i][1],index[i][2],index[i][3],index[i][4],index[i][5], useCharge, fasta));
		}	
	}
	
	private  AnnotatedSpectrum createDummy(){
		AnnotatedSpectrum dummy = new AnnotatedSpectrum();
		dummy.peptide = "DUMMY";
		dummy.protein = "NO-PROTEIN";
		dummy.score = 0;
		dummy.getAnnotation().put("fdr", 1.0);
		dummy.getAnnotation().put("pepfdr", 1);
		dummy.getAnnotation().put("specCount", 0);
		return dummy;
	}
	
	public void getMultiRunPeptideStat(){
		int total = this.results.size();
		Set<String> combinedPeps = new HashSet<String>();
		for(int i = 0; i < results.size(); i++){
			TDAStat result = results.get(i);
			combinedPeps.addAll(result.peptideMap.keySet());
		}
		Map<String, List<AnnotatedSpectrum>> pepMap = new HashMap<String, List<AnnotatedSpectrum>>(combinedPeps.size());

		for(Iterator<String> it = combinedPeps.iterator(); it.hasNext();){
			List<AnnotatedSpectrum> results = new ArrayList<AnnotatedSpectrum>(total);
			pepMap.put(it.next(), results);
		}
		AnnotatedSpectrum dummy = createDummy();
		for(Iterator<String> it = combinedPeps.iterator(); it.hasNext();){
			String pep = it.next();
			for(int i = 0; i < results.size(); i++){
				TDAStat result = results.get(i);
				List<AnnotatedSpectrum> results = pepMap.get(pep);
				if(result.peptideMap.containsKey(pep)){
					results.add(result.peptideMap.get(pep));
				}else{
					results.add(dummy);
				}
				
			}
		}
		
		for(Iterator<String> it = pepMap.keySet().iterator(); it.hasNext();){
			String pep = it.next();
			String prot = "";
			Collection values = pepMap.get(pep);
			Iterator it2 = values.iterator();
			System.out.print(pep);
			int count=1;
			while(it2.hasNext()){
				AnnotatedSpectrum rest = (AnnotatedSpectrum)it2.next();
				if(!rest.protein.equals("NO-PROTEIN")){
					prot = rest.protein;
				}
				System.out.print("\t" + rest.score + "\t" + rest.getAnnotation().get("specCount"));
				count++;
			}
			System.out.print("\t" + prot);
			System.out.println();
			
		}
		
	}
		
	
	public void getPairWiseSim(){
		TDAStat results1 = this.results.get(0);
		TDAStat results2 = this.results.get(1);
		Map<String, AnnotatedSpectrum> peps1 = results1.peptideMap;
		Map<String, AnnotatedSpectrum> peps2 = results2.peptideMap;
		MZXMLReader reader1 = readers.get(0);
		MZXMLReader reader2 = readers.get(1);
		for(Iterator<String> it = peps1.keySet().iterator(); it.hasNext();){
			String pep = it.next();
			if(peps2.containsKey(pep)){
				Spectrum libSpect = this.lib.getSpectra(pep).get(0);
				libSpect.filterPeaks(30);
				libSpect.sqrtSpectrum();
				AnnotatedSpectrum psm1 = peps1.get(pep);
				AnnotatedSpectrum psm2 = peps2.get(pep);
				Spectrum swath1 = reader1.getSpectrum(psm1.scanNumber);
				swath1.windowFilterPeaks2(15, 25);
				swath1.mergePeaks(swath1, 0.05);
				swath1.sqrtSpectrum();
				Spectrum swath2 = reader2.getSpectrum(psm2.scanNumber);
				swath2.mergePeaks(swath2, 0.05);
				swath2.windowFilterPeaks2(15, 25);
				swath2.sqrtSpectrum();
				Spectrum project1 = swath1.project(libSpect, 0.05);
				Spectrum project2 = swath2.project(libSpect, 0.05);
				System.out.println(pep + "\t" + psm1.protein + "\tsim:\t" + psm1.score + "\t" + psm2.score 
						+"\t" + libSpect.projectedCosine(swath1, 0.05) + "\t" + libSpect.projectedCosine(swath2, 0.05)
						+ "\t" + swath1.cosine(swath2, 0.05) + "\t" + libSpect.getPeak().size()
						+ "\t" + project1.cosine(project2, 0.05));
			}
		}
	}
	
	public void getBackGroundPairSim(){
		MZXMLReader reader1 = readers.get(0);
		MZXMLReader reader2 = readers.get(1);
		int maxCount = 10000000;
		int size = reader1.getSpectrumCount();
		int size2 = reader2.getSpectrumCount();
		for(int i = 0; i < 10000; i++){
			Spectrum s = reader1.getSpectrum((int)(Math.random()*size+1));
			s.windowFilterPeaks2(15, 25);
			s.mergePeaks(s, 0.05);
			s.sqrtSpectrum();
			for(int j = 0; j < 500; j++){
				Spectrum s2 = reader2.getSpectrum((int)(Math.random()*size2+1));
				if(Math.abs(s.parentMass - s2.parentMass) < 2){
					s2.windowFilterPeaks2(15, 25);
					s2.mergePeaks(s, 0.05);
					s2.sqrtSpectrum();
					System.out.println(s.scanNumber + "\t" + s2.scanNumber + "\tsim:\t" + s.cosine(s2, 0.05));
				}
			}
		}
		
	}
	
	
	public static void testPairSim(){
		String swath1 = "../mixture_linked/msdata/UPS12_Human/18470_REP3_40fmol_UPS1_500ng_HumanLysate_SWATH_1.mzXML";
		String swath2 = "../mixture_linked/msdata/UPS12_Human/18474_REP3_40fmol_UPS2_500ng_HumanLysate_SWATH_1.mzXML";
		String result1 = "../mixture_linked/SWATH/Swath_Human_searches/Swathmsplit/UPS_Human_REP3withLysatelib_msplit_1_pep.txt";
		String result2 = "../mixture_linked/SWATH/Swath_Human_searches/Swathmsplit/UPS_Human_REP3withLysatelib_msplit_3_pep.txt";
		String libFile = "../mixture_linked/UPS_Human_lysateREP123_IDA_1pepFDR_plusDecoy2.mgf";
		MultiSWATHAnalyzer multiSwath = new MultiSWATHAnalyzer(swath1, swath2, result1, result2, libFile);
		//multiSwath.getPairWiseSim();
		multiSwath.getBackGroundPairSim();
		
	}
	
	public static void testMultiSWATH(String[] files){
//		String fasta = "../mixture_linked/database/Human_allproteins_plusDecoy.fasta";
//		String result0  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_0.txt";
//		String result1  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_1.txt";
//		String result2  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_2.txt";
//		String result3  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_3_pep.txt";
//		String result4  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_4_pep.txt";
//		String result5  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_5_pep.txt";
//		String result6  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_6_pep.txt";
//		String result7  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_7_pep.txt";
//		String result8  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/APSWATHSwathmsplit/APSWATH_withREP3combinedLib_msplit_8_pep.txt";

		String fasta = "../mixture_linked/database/Human_allproteins_plusDecoy.fasta";
//		String result0  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/PeakView/APLib/Swath_EIF4aJune7-Biorep1 - Peptides.txt";
//		String result1  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/PeakView/APLib/Swath_EIF4aJune7-Biorep2 - Peptides.txt";
//		String result2  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/PeakView/APLib/Swath_EIF4aJune7-Biorep3 - Proteins - Peptides.txt";
		String result3  = "../mixture_linked/Swath_MEPCE_Biorep1_withLysatelib_withRTcalc_5minWidth_quant - Peptides.txt";
		String result4  = "../mixture_linked/Swath_GFP_Biorep1_withLysatelib_withRTcalc_5minWidth_quant - Peptides.txt";
		String result5  = "../mixture_linked/Swath_EIF4A2_Biorep3_withLysatelib_withRTcalc_5minWidth_quant - Peptides.txt";
//		String result6  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/PeakView/APLib/Swath_MEPCEJune7-Biorep1 - Peptides.txt";
//		String result7  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/PeakView/APLib/Swath_MEPCEJune7-Biorep2 - Peptides.txt";
//		String result8  = "../mixture_linked/SWATH/APSWATH/PeptideIDs/PeakView/APLib/Swath_MEPCEJune7-Biorep3 - Peptides.txt";
		
		int[] swathIndex = new int[]{1, 4, 6, 8,32,-1};
		int[] msgfdbInd = new int[]{1, 7, 6, 8,11,1};
		int[] peakvwInd = new int[]{-1, 1, 3, 0, 5,1};
//<--------------------------- UPS samples -------------------------------->		
//		String fasta = "../mixture_linked/database/UPS2.fasta";
//		String result0  = "../mixture_linked/SWATH/UPS_only/Swathmsplit/UPS_only_UPSEcoliREP234lib_msplit_2.txt";
//		String result1  = "../mixture_linked/SWATH/UPS_only/Swathmsplit/UPS_only_UPSEcoliREP234lib_msplit_3.txt";
//		String result2  = "../mixture_linked/SWATH/UPS_only/Swathmsplit/UPS_only_UPSEcoliREP234lib_msplit_0.txt";
//		String result3  = "../mixture_linked/SWATH/UPS_only/Swathmsplit/UPS_only_UPSEcoliREP234lib_msplit_1.txt";
//		String result4  = "../mixture_linked/SWATH/UPS_only/IDA/UPS_REP3_IDA_msgfdb_2.txt";
//		String result5  = "../mixture_linked/SWATH/UPS_only/IDA/UPS_REP3_IDA_msgfdb_3.txt";
//		String result6  = "../mixture_linked/SWATH/UPS_only/IDA/UPS_REP3_IDA_msgfdb_0.txt";
//		String result7  = "../mixture_linked/SWATH/UPS_only/IDA/UPS_REP3_IDA_msgfdb_1.txt";
//		String result8  = "../mixture_linked/SWATH//UPS_only/PeakView/CombinedDDALib/18511_REP3_40fmol_UPS1_PeakViewExtractDDAlib - Peptides.txt";
//		String result9  = "../mixture_linked/SWATH/UPS_only/PeakView/CombinedDDALib/1852_REP3_400fmol_UPS1_PeakViewExtractDDAlib - Peptides.txt";
//		String result10  = "../mixture_linked/SWATH/UPS_only/PeakView/CombinedDDALib/18507_REP3_40fmol_UPS2_PeakViewExtractDDLib - Peptides.txt";
//		String result11 = "../mixture_linked/SWATH/UPS_only/PeakView/CombinedDDALib/18509_REP3_400fmol_UPS2_PeakViewExtractDDAlib - Peptides.txt";

//<--------------------------- UPS Ecoli samples -------------------------------->		
//		String fasta = "../mixture_linked/database/UPS_plusEcoli_plusDecoy.fasta";
//		String result0  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPSEcoli_libREp234_msplit_4.txt";
//		String result1  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPSEcoli_libREp234_msplit_5.txt";
//		String result2  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPSEcoli_libREp234_msplit_7.txt";
//		String result3  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/SWATH/UPSEcoli_libREp234_msplit_8.txt";
//		String result4  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/IDA/UPSEcoli_REP3_newStock2_msgfdb_1.txt";
//		String result5  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/IDA/UPSEcoli_REP3_newStock2_msgfdb_2.txt";
//		String result6  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/IDA/UPSEcoli_REP3_newStock2_msgfdb_4.txt";
//		String result7  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/IDA/UPSEcoli_REP3_newStock2_msgfdb_5.txt";
//		String result8  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/PeakView/CombinedDDALib/18484_REP3_1ug_Ecoli_PeakViewExtractDDAlib - Peptides.txt";
//		String result9  = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/PeakView/CombinedDDALib/18488_REP3_40fmol_UPS1_1ug_Ecoli_PeakViewExtractDDAlib - Peptides.txt";
//		String result10 = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/PeakView/CombinedDDALib/18492_REP3_40fmol_UPS2_1ugEcoli_PeakViewExtractDDALib - Peptides.txt";
//		String result11 = "../mixture_linked/SWATH/UPS_EcoliREP2_REP3_searches/PeakView/CombinedDDALib/18494_REP3_400fmol_UPS2_1ug_Ecoli_PeakViewExtractDDAlib - Peptides.txt";
		
//<-------------------UPS  + Human samples--------------------------->	
//		String fasta = "../mixture_linked/database/Human_allproteins_plusDecoy.fasta";
//		String result0  = "../mixture_linked/SWATH/Swath_Human_searches/Swathmsplit/UPS_Human_REP3withLysatelib_msplit_0.txt";
//		String result1  = "../mixture_linked/SWATH/Swath_Human_searches/Swathmsplit/UPS_Human_REP3withLysatelib_msplit_1.txt";
//		String result2  = "../mixture_linked/SWATH/Swath_Human_searches/Swathmsplit/UPS_Human_REP3withLysatelib_msplit_3.txt";
//		String result3  = "../mixture_linked/SWATH/Swath_Human_searches/Swathmsplit/UPS_Human_REP3withLysatelib_msplit_4.txt";
//		String result4  = "../mixture_linked/SWATH/Swath_Human_searches/IDAMsgfdb/UPS_Human_IDAREP3_msgfdb_0.txt";
//		String result5  = "../mixture_linked/SWATH/Swath_Human_searches/IDAMsgfdb/UPS_Human_IDAREP3_msgfdb_1.txt";
//		String result6  = "../mixture_linked/SWATH/Swath_Human_searches/IDAMsgfdb/UPS_Human_IDAREP3_msgfdb_3.txt";
//		String result7  = "../mixture_linked/SWATH/Swath_Human_searches/IDAMsgfdb/UPS_Human_IDAREP3_msgfdb_4.txt";
//		String result8  = "../mixture_linked/SWATH/Swath_Human_searches/PeakView/CombinedDDALib/18468_REP3_500ng_HumanLysate_PeakViewExtractDDAlib - Peptides.txt";
//		String result9  = "../mixture_linked/SWATH/Swath_Human_searches/PeakView/CombinedDDALib/18470_REP3_40fmol_UPS1_500ng_HumanPeakViewExtractDDAlib - Peptides.txt";
//		String result10  = "../mixture_linked/SWATH/Swath_Human_searches/PeakView/CombinedDDALib/18474_REP3_40fmol_UPS2_500ng_HumanLysate_PeakViewExtractDDAlib - Peptides.txt";
//		String result11 = "../mixture_linked/SWATH/Swath_Human_searches/PeakView/CombinedDDALib/18476_REp3_400fmol_UPS2_500ng_Human_PeakViewExtractDDAlib - Peptides.txt";
//		files = new String[]{result0, result1, result2, result3, result4, result5, result6, result7, result8, result9, result10, result11};
//		int[][] index = new int[][]{swathIndex, swathIndex, swathIndex, swathIndex, msgfdbInd, msgfdbInd, msgfdbInd, msgfdbInd, peakvwInd, peakvwInd, peakvwInd, peakvwInd};
		files = new String[]{result3, result4, result5};
		int[][] index = new int[][]{peakvwInd, peakvwInd, peakvwInd};
//		files = new String[]{result0, result1, result2, result3};
//		int[][] index = new int[][]{swathIndex, swathIndex, swathIndex, swathIndex};
		MultiSWATHAnalyzer multiSwath = new MultiSWATHAnalyzer(files, index, fasta);
		multiSwath.getMultiRunPeptideStat();
	}
	
	
	public static void testMultiSWATHDir(){
		String resultDir = "..//mixture_linked//Emily_toni_Exosomes/DIA_result/";
		File dir = new File(resultDir);
		String[] files = dir.list();
		//int[] msgfdbInd =  new int[]{1, 7, 6, 8,11,1};
		int[] msplitDIAInd = new int[]{1, 4, 6, 8,32,-1};
		int[][] index = new int[files.length][];
		for(int i = 0; i < index.length; i++){
			index[i] = msplitDIAInd;
		}
		for(int i = 0; i < index.length; i++){
			files[i] = resultDir + files[i];
		}
		String fasta = "../mixture_linked/Emily_toni_Exosomes/mergedSequence/1798464a5a4549d68de1c131e8500d1f.fasta";
		MultiSWATHAnalyzer multiSwath = new MultiSWATHAnalyzer(files, index, fasta);
		multiSwath.getMultiRunPeptideStat();
	}
	
	public static void main(String[] args){
		//testPairSim();
		testMultiSWATH(args);
		//testMultiSWATHDir();
	}
	

}
