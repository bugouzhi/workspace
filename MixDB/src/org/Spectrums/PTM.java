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
	public double ptmMass = 0.0;
	private int[] targetedPos = new int[0]; //targeted position
	private char[] targetedRes = new char[0];
	
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
			if(pos < peptide.length() && pos != CrossLinker.ANYPOSITION && pos != CrossLinker.CTERM &&
					pos != CrossLinker.NTERM &&
					(targetedRes[i] == peptide.charAt(pos) 
							|| targetedRes[i] == CrossLinker.ANYRESIDUE)){
				positions.add(pos+1);
			}
			if(pos == CrossLinker.CTERM){ 
					//&& peptide.charAt(peptide.length()-1) == targetedRes[i]){
				positions.add(peptide.length());
			}
			if(pos == CrossLinker.NTERM){
				positions.add(1);
			}
			if(pos == CrossLinker.ANYPOSITION){
				for(int j = 0; j < peptide.length(); j++){
					if(peptide.charAt(j) == targetedRes[i] ||
							targetedRes[i] == CrossLinker.ANYPOSITION){
						positions.add(j+1);
					}
				}
			}
		}
		List<Integer> ret = new ArrayList<Integer>();
		ret.addAll(positions);
		//System.out.println(ret.size());
		return ret;
	}
	
	public String toString(){
		return "{" + this.ptmMass + "\t" 
				+ Arrays.toString(this.targetedPos) 
				+ "\t" + Arrays.toString(this.targetedRes);
	}
	
	public static List<PTM[]> generatePTMList(List<PTM> ptms, int maxPTM){
		List<PTM[]> ptmList = new ArrayList<PTM[]>();
		if(maxPTM > 1){
			List<PTM[]> ptmSubList = generatePTMList(ptms, maxPTM-1);
			for(int i = 0; i < ptmSubList.size(); i++){
				PTM[] ptmCom = ptmSubList.get(i);
				if(ptmCom.length == maxPTM - 1){
					insertOnePTM(ptmCom, ptms, ptmSubList);
				}
			}
			//System.out.println("ptm size: " + ptmSubList.size());
			return ptmSubList;
		}else{
			for(Iterator<PTM> it = ptms.iterator(); it.hasNext();){
				ptmList.add(new PTM[]{it.next()});
			}
			return ptmList;
		}
	}
	
	public static List<PTM[]> insertOnePTM(PTM[] ptm, List<PTM> ptms, List<PTM[]> ptmList){
		for(int i = 0; i < ptms.size(); i++){
			PTM[] newPtms = Arrays.copyOf(ptm, ptm.length+1);
			newPtms[newPtms.length-1] = ptms.get(i);
			ptmList.add(newPtms);
		}
		return ptmList;
	}
}
	
	