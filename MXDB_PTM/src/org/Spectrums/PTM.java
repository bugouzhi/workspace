package org.Spectrums;

import java.util.ArrayList;
import java.util.HashSet;
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
		return ret;
	}
	
	
}
