package org.Spectrums;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A cross-linker is a chemical agent that covalentlly link two amino acid
 * residues that are spatially close together
 * @author Jian Wang
 *
 */
public class CrossLinker {
	public static int ANYPOSITION = -1;
	public static int NTERM = 0;
	public static int CTERM = -2;
	public static char ANYRESIDUE = '*';
	private static int MONOFUNCTIONAL = 1;
	private static int BIFUNCTIONAL = 2;
	private int type; //type of cross-linker, bi-functional or mono-functional so-far, can add more later
	private double[] linkerMass = new double[]{0.0}; //mass of the cross-linker
	private double[] linkerMassOffSets = new double[]{0.0}; //mass offset added to linked peptides
	private int[] targetedPos1 = new int[0]; //targeted position
	private int[] targetedPos2 = new int[0];
	private char[] targetedRes1 = new char[0];
	private char[] targetedRes2 = new char[0];
	private int[] targetPos=new int[0];
	private char[] targetRes=new char[0];
	
	public static CrossLinker DSP = new CrossLinker(Mass.DSPLINKER_MASS, new char[]{'K'});
	public static CrossLinker DSS = new CrossLinker(Mass.DSSLINKER_MASS, new char[]{'K'});

	//homobifunctional crosslinkers
	public CrossLinker(double linkerMass, char[] targetedPositions){
		this.linkerMassOffSets = new double[]{linkerMass};
		this.targetedPos1 = new int[]{CrossLinker.ANYPOSITION};
		this.targetedPos2 = new int[]{CrossLinker.ANYPOSITION};		
		this.targetedRes1 = targetedPositions;
		this.targetedRes2 = targetedPositions;
	}
	
	
	public CrossLinker(double[] linkerMasses, char[] targetedPositions){
		this.linkerMassOffSets = linkerMasses;
		this.targetedPos1 = new int[]{CrossLinker.ANYPOSITION};
		this.targetedPos2 = new int[]{CrossLinker.ANYPOSITION};		
		this.targetedRes1 = targetedPositions;
		this.targetedRes2 = targetedPositions;
	}
	
	public CrossLinker(double linkerMass, int[] targetPos1, int[] targetPos2, 
			char[] targetRes1, char[] targetRes2){
		this.linkerMassOffSets = new double[]{linkerMass};
		this.targetedPos1 = targetPos1;
		this.targetedPos2 = targetPos2;
		this.targetedRes1 = targetRes1;
		this.targetedRes2 = targetRes2;
	}
	
	public double getLinkerMassOffSet(){
		return this.linkerMassOffSets[0];
	}
	
	public double[] getLinkerMassOffSets(){
		return this.linkerMassOffSets;
	}
	
	
	//to be implemented
	//the mass of linker only linked one-side
	//to-be-implemented
	public double getLinkerDangleMass(){
		return 0.0;
	}
	
	private void mergeLinkerSpecificity(){
		List<Integer> list1 = new ArrayList();
		List<Character> list2 = new ArrayList();
		for(int i = 0; i < this.targetedPos1.length; i++){
			int position1 = this.targetedPos1[i];	
			for(int j = 0; j < this.targetedPos2.length; j++){
				int position2 = this.targetedPos2[j];	
			}
		}
	}
	//get possible linking position for peptide 1
	public List<Integer> getLinkerPositions1(String peptide){
		return getLinkerPositions(peptide, this.targetedPos1, this.targetedRes1);
	}
	
	//get possible linking position for peptide 2
	public List<Integer> getLinkerPositions2(String peptide){
		return getLinkerPositions(peptide, this.targetedPos2, this.targetedRes2);
	}
	
	//get possible linking position for peptide regardless whehter it is one or two
	public List<Integer> getLinkerPositions(String peptide){
		return null;
	}
	
	public List<Integer> getLinkerPositions(String peptide, int[] targetedPos, char[] targetedRes){
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
	
	
	public List<LinkedPeptide> crossLinkPeptides(Peptide p1, Peptide p2, int minCharge, int maxCharge){
		List<Integer> positions = getLinkerPositions1(p1.getPeptide());
		List<Integer> positions2 = getLinkerPositions2(p1.getPeptide());
		List<LinkedPeptide> toBeAdded = new ArrayList<LinkedPeptide>();
		for(int c = minCharge; c <= maxCharge; c++){
			for(int m = 0; m < positions.size(); m++){
				for(int n = 0; n < positions2.size(); n++){
					LinkedPeptide linked = new LinkedPeptide(p1, p2, c, positions.get(m)+1, positions2.get(n)+1);
					toBeAdded.add(linked);
				}
			}
		}
		return toBeAdded;
	}
	
	public boolean isLinkedPosition(char aa){
		for(int i = 0; i < targetedRes1.length; i++){
			if(targetedRes1[i] == aa){
				return true;
			}
		}
		return false;
	}
}
