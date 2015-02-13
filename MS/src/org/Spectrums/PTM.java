package org.Spectrums;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * A object represent a specific PTM on a peptide
 * @author jian wang
 *
 */
public class PTM {
	public static int ANYPOSITION = -1;
	public static int NTERM = 0;
	public static int CTERM = -2;
	public static char ANYRESIDUE = '*';
	private double ptmMass = 0.0;
	private int[] targetedPos = new int[0]; //targeted position
	private char[] targetedRes = new char[0];
	private String name;  //name of the PTM
	
	public PTM(double ptmMass, int[] targetPos, char[] targetRes){
		this.ptmMass = ptmMass;
		this.targetedPos = targetPos;
		this.targetedRes = targetRes;
	}
	
	public List<Integer> getPTMPositions(String peptide){
		return getPTMPositions(peptide, this.targetedPos, this.targetedRes);
	}
	
	private List<Integer> getPTMPositions(String peptide, int[] targetedPos, char[] targetedRes){
		Set<Integer> positions = new HashSet<Integer>();
		for(int i = 0; i < targetedPos.length; i++){
			int pos = targetedPos[i];
			if(pos < peptide.length() && pos != CrossLinker.ANYPOSITION && //pos != CrossLinker.CTERM &&  
					(targetedRes[i] == peptide.charAt(pos) 
							|| targetedRes[i] == CrossLinker.ANYRESIDUE)){
				positions.add(pos);
			}
			if(pos == CrossLinker.CTERM 
					&& peptide.charAt(peptide.length()-1) == targetedRes[i]){
				positions.add(peptide.length()-1);
			}
			if(pos == CrossLinker.ANYPOSITION){
				for(int j = 0; j < peptide.length(); j++){
					if(peptide.charAt(j) == targetedRes[i] ||
							targetedRes[i] == CrossLinker.ANYPOSITION){
						positions.add(j);
					}
				}
			}
		}
		List<Integer> ret = new ArrayList<Integer>();
		ret.addAll(positions);
		//System.out.println("ret-size: " + ret.size());
		return ret;
	}

	public double getPtmMass() {
		return ptmMass;
	}

	public void setPtmMass(double ptmMass) {
		this.ptmMass = ptmMass;
	}
	
	/**
	 * Parse PTM files to generate the possible ptms allowed for this search
	 * @param precursorsFile
	 * @return
	 */
	public static List<PTM[]> parsePTMs(String ptmFile){
		List<PTM[]> ptmList = new ArrayList();
		List<String> lines = Utils.FileIOUtils.createListFromFile(ptmFile);
		int maxPtm = 0;
		List<PTM>ptms = new ArrayList<PTM>();
//		this.ptms.add(new PTM(28.03, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(32.06, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(4.01, new int[]{PTM.CTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(16.01, new int[]{PTM.ANYPOSITION}, new char[]{'M'}));
		for(int i = 0; i < lines.size(); i++){
			String[] line = lines.get(i).split(",");
			if(line[0].startsWith("#")){
				continue;
			}else if(line[0].equals("maxPTM")){
				maxPtm = Integer.parseInt(line[1]);
			}else{
				PTM ptm = PTM.parsePTM(line);
				ptms.add(ptm);
			}
		}
		ptmList = PTM.generatePTMList(ptms, 1);	
		System.out.println("size: " + ptmList.size());
		//ptmList.add(new PTM[]{});
		System.out.println("Finish getting ptms info, ptm size: " + ptmList.size());
		return ptmList;
	}
	/**
	 * Parse PTMs from XML files which is generated from the param.xml in ProteoSAfe
	 * @param xmlFile
	 * @return
	 */
	public static List<PTM[]> parsePTMFromXML(String xmlFile){
		List<PTM[]> ptmList = new ArrayList();
		List<String> lines = Utils.FileIOUtils.createListFromFile(xmlFile);
		int maxPtm = 0;
		List<PTM>ptms = new ArrayList<PTM>();
//		this.ptms.add(new PTM(28.03, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(32.06, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(4.01, new int[]{PTM.CTERM}, new char[]{PTM.ANYRESIDUE}));
//		this.ptms.add(new PTM(16.01, new int[]{PTM.ANYPOSITION}, new char[]{'M'}));
		for(int i = 0; i < lines.size(); i++){
			String line = lines.get(i);
			if(line.startsWith("<parameter name=\"ptm.mods\">")){
				maxPtm = Integer.parseInt(line.split("[<>]")[2]);
			}else if(line.startsWith("<parameter name=\"ptm.custom_PTM\">")){
				String[] tokens = line.split("[<>]");
				//System.out.println("tokens: 1" + tokens[2]);
				String[] tokens2 = tokens[2].split(",");
				String[] ptmStr = new String[4];
				ptmStr[0] = tokens2[0];
				ptmStr[1] = tokens2[1];
				ptmStr[2]="any";
				if(tokens2[2].contains("nterm")){
					ptmStr[2] = "Nterm";
				}
				if(tokens2[2].contains("cterm")){
					ptmStr[2] = "Cterm";
				}
				ptmStr[3]="opt";
				PTM ptm = PTM.parsePTM(ptmStr);
				ptms.add(ptm);
			}
		}
		ptmList = PTM.generatePTMList(ptms, 1);	
		System.out.println("size: " + ptmList.size());
		//ptmList.add(new PTM[]{});
		System.out.println("Finish getting ptms info, ptm size: " + ptmList.size());
		return ptmList;
	}
	
	private static PTM parsePTM(String[] tokens){
		double PtmMass = Double.parseDouble(tokens[0]);
		String residues = tokens[1];
		char[] targetedRes = new char[residues.length()];
		int[] targetedPos = new int[residues.length()];
		int pos=0;
		for(int i = 0; i < residues.length(); i++){
			if(residues.charAt(i) == '*'){
				targetedRes[i] = PTM.ANYRESIDUE;
			}else{
				targetedRes[i] = residues.charAt(i);
			}
		}
		
		if(tokens[2].equals("N-term")){
			pos = PTM.NTERM;
		}else if(tokens[2].equals("C-term")){
			pos = PTM.CTERM;
		}else if(tokens[2].equals("any")){
			pos = PTM.ANYPOSITION;
		}else{
			pos = Integer.parseInt(tokens[2]);
		}
		targetedPos = new int[]{pos};
		System.out.println("parsing ptm");
		return new PTM(PtmMass, targetedPos, targetedRes);
	}
	
	public static List<PTM[]> generatePTMList(List<PTM> ptms, int maxPTM){
		List<PTM[]> ptmList = new ArrayList<PTM[]>();
		if(maxPTM > 1){
			List<PTM[]> ptmSubList = generatePTMList(ptms, maxPTM-1);
			for(int i = 0; i < ptmSubList.size(); i++){
				PTM[] ptmCom = ptmSubList.get(i);
				if(ptmCom.length == maxPTM - 1){
					insertOnePTM(ptmCom, ptms, ptmSubList, i % ptms.size());
					System.out.println("beginInd " + i % ptms.size());
				}
			}
			System.out.println("ptm size: " + ptmSubList.size() + "\t" + maxPTM);
			return ptmSubList;
		}else{
			for(Iterator<PTM> it = ptms.iterator(); it.hasNext();){
				ptmList.add(new PTM[]{it.next()});
			}
			System.out.println("ptm size: " + ptmList.size() + "\t" + maxPTM);
			return ptmList;
		}
	}
	
	
	
	public static List<PTM[]> insertOnePTM(PTM[] ptm, List<PTM> ptms, List<PTM[]> ptmList, int beginInd){
		for(int i = beginInd; i < ptms.size(); i++){
			PTM[] newPtms = Arrays.copyOf(ptm, ptm.length+1);
			newPtms[newPtms.length-1] = ptms.get(i);
			ptmList.add(newPtms);
		}
		return ptmList;
	}
	
	
}
