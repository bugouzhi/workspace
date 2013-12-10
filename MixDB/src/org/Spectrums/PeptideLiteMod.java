package org.Spectrums;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Represents peptide with PTMs
 * @author Jian
 *
 */
public class PeptideLiteMod extends PeptideLite{
	private static final long serialVersionUID = 1231234L;
	private int[] ptmPos;
	private double[] ptmMasses;
	public PeptideLiteMod(int beginInd, int endInd) {
		super(beginInd, endInd);
		this.ptmPos = new int[]{};
		this.ptmMasses = new double[]{};
	}
	
	//insert a ptm to the peptide
	//if an residue is already modified, will not add a second
	//modification to it
	public PeptideLiteMod insertMod(int pos, double mass) {
		for(int i = 0; i < this.ptmPos.length; i++){
			if(ptmPos[i] == pos){
				return null;
			}
		}
		//position orders are enforced to eliminate cases where
		//multi combo of same mods set added in different order
		//are not consider during the addition of modifications
		for(int i = this.ptmMasses.length-1; i >= 0; i--){
			//System.out.println(mass + "\t" + ptmMasses[i]);
			if(pos <= ptmPos[i]){
				return null;
			}
		}
		PeptideLiteMod modPep = super.insertMod(this.getBeginInd(), this.getEndInd());
		modPep.ptmPos = Arrays.copyOf(this.ptmPos, this.ptmPos.length+1);
		modPep.ptmMasses = Arrays.copyOf(this.getPtmMasses(), this.getPtmMasses().length+1);
		modPep.ptmPos[modPep.getPtmPos().length-1]=pos;
		modPep.getPtmMasses()[modPep.getPtmMasses().length-1]=mass;
		return modPep;
	}

	public int[] getPtmPos() {
		return ptmPos;
	}
	public void setPtmPos(int[] ptmPos) {
		this.ptmPos = ptmPos;
	}
	public double[] getPtmMasses() {
		return ptmMasses;
	}
	public void setPtmMasses(double[] ptmMasses) {
		this.ptmMasses = ptmMasses;
	}
	
	public String toString(){
		return getModPep();
	}
	
	public String getModPep(){
		StringBuffer buff = new StringBuffer(super.toString());
		int offset=0;
		for(int i = 0; i < this.ptmPos.length; i++){
			String ptm;
			if(this.ptmMasses[i]>0){
				ptm = String.format("+%.3f", this.ptmMasses[i]);
			}else{
				ptm = String.format("%.3f", this.ptmMasses[i]);
			}
			buff.insert(this.ptmPos[i]+offset, ptm);
			offset+=ptm.length();
		}
		return buff.toString();
	}
	
	//add PTMs to peptides, modified peptides can have up to maxPTM
	public static List<PeptideLite> insertPTM(List<PeptideLite> pepList, List<PTM> ptms, int maxPTM){
		List<PeptideLite> modPeps = new ArrayList<PeptideLite>();
		if(maxPTM > 1){
			List<PeptideLite> modPepsSub = insertPTM(pepList, ptms, maxPTM-1);
			modPeps.addAll(modPepsSub);
			List<PeptideLite> modPepsLast = insertOnePTM(modPepsSub, ptms, maxPTM);
			modPeps.addAll(modPepsLast);
			return modPeps;
		}else{
			return insertOnePTM(pepList, ptms, maxPTM);
		}
	}
	
	//add exactly one PTM from the ptm list to peptide list, mod peptide can have up to maxPTM, peptides already with
	//these PTMs will not be added any further PTMs
	//return new modified peptides
	public static List<PeptideLite> insertOnePTM(List<PeptideLite> pepList, List<PTM> ptms, int maxPTM){
		List<PeptideLite> modPeps = new ArrayList<PeptideLite>();
		for(int i = 0; i < ptms.size(); i++){
			List<PeptideLite> modPepOne = insertOnePTM(pepList, ptms.get(i), maxPTM);
			modPeps.addAll(modPepOne);
		}
		return modPeps;
	}
	
	//add ptm to the peptide list
	public static List<PeptideLite> insertOnePTM(List<PeptideLite> pepList, PTM ptm, int maxPTM){
		List<PeptideLite> modPeps = new ArrayList<PeptideLite>();
		for(int i = 0; i < pepList.size(); i++){
			PeptideLite pep = pepList.get(i);
			List<Integer> ptmpos = ptm.getPTMPositions(pep.getPep());
			for(int j = 0; j < ptmpos.size(); j++){
				int pos = ptmpos.get(j);
				PeptideLiteMod mod = pep.insertMod(ptmpos.get(j), ptm.ptmMass);
				if(mod != null){
					modPeps.add(mod);
				}
			}
		}
		return modPeps;
	}
	
	public static void testinsertPTM(){
		String testProtein = "KYAFSLTTFSPNGKLVQIEYALNAVNAGVTSVGIKATDGVVLATEKKPTSELAIGAS";
		List<PeptideLite> pepList = new ArrayList<PeptideLite>();
		pepList.add(new PeptideLite(0, 12, testProtein, 2));
		pepList.add(new PeptideLite(5, 20, testProtein, 3));
		pepList.add(new PeptideLite(15, 28, testProtein, 2));
		pepList.add(new PeptideLite(20, 35, testProtein, 3));
		pepList.add(new PeptideLite(36, 45, testProtein, 2));
		for(int i = 0; i < pepList.size(); i++){
			System.out.println("peps: " + pepList.get(i));
		}
		List<PTM> ptms = new ArrayList<PTM>();
		ptms.add(new PTM(16.01, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
		ptms.add(new PTM(40.01, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
		ptms.add(new PTM(24.01, new int[]{PTM.CTERM}, new char[]{PTM.ANYRESIDUE}));
		List<PeptideLite> modPeps = insertPTM(pepList, ptms, 2);
		System.out.println("modPeps has size: " + modPeps.size());
		for(int i = 0; i < modPeps.size(); i++){
			System.out.println("peps: " + modPeps.get(i));
		}

	}
	
	public static void testGeneratPTMs(){
		List<PTM> ptms = new ArrayList<PTM>();
		ptms.add(new PTM(16.01, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
		ptms.add(new PTM(40.01, new int[]{PTM.NTERM}, new char[]{PTM.ANYRESIDUE}));
		ptms.add(new PTM(24.01, new int[]{PTM.CTERM}, new char[]{PTM.ANYRESIDUE}));
		List<PTM[]> ptmList = PTM.generatePTMList(ptms, 2);
		for(int i = 0; i < ptmList.size(); i++){
			System.out.println("Mods are: " + Arrays.toString(ptmList.get(i)));
		}
	}
	
	public static void main(String[] args){
		//testinsertPTM();
		testGeneratPTMs();
	}
		
}
