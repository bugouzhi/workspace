package org.Spectrums;
import java.io.Serializable;

import sequences.FastaSequence;

/**
 * A simple memory efficient implementation of peptides
 * @author Jian Wang
 *
 */
public class PeptideLite implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 230483204832L;
	private static String EMPTYSTR = "";
	private int beginInd;
	private int endInd;
	private int charge;
	private String protein;
	private FastaSequence fastaseq;
	
	public String getProtein() {
		return protein;
	}

	public void setProtein(String protein) {
		this.protein = protein;
	}
	
	public PeptideLite(int beginInd, int endInd, String proteins, int charge){
		this.beginInd= beginInd;
		this.endInd = endInd;
		this.protein = proteins;
		this.charge = charge;
	}
	
	public PeptideLite(int beginInd, int endInd, FastaSequence seq){
		this(beginInd, endInd, EMPTYSTR, 1);
		this.fastaseq = seq;
	}
	
	//created a modify peptides from this peptide
	public PeptideLiteMod insertMod(int pos, double mass) {
		PeptideLiteMod modPep = new PeptideLiteMod(this.getBeginInd(), this.getEndInd());
		modPep.setFastaseq(this.fastaseq);
		modPep.setProtein(this.protein);
		modPep.setPtmPos(new int[]{pos});
		modPep.setPtmMasses(new double[]{mass});
		return modPep;
	}
	
	public static String getEMPTYSTR() {
		return EMPTYSTR;
	}

	public static void setEMPTYSTR(String emptystr) {
		EMPTYSTR = emptystr;
	}

	public int getBeginInd() {
		return beginInd;
	}

	public void setBeginInd(int beginInd) {
		this.beginInd = beginInd;
	}

	public int getEndInd() {
		return endInd;
	}

	public void setEndInd(int endInd) {
		this.endInd = endInd;
	}

	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public FastaSequence getFastaseq() {
		return fastaseq;
	}

	public void setFastaseq(FastaSequence fastaseq) {
		this.fastaseq = fastaseq;
	}

	public PeptideLite(int beginInd, int endInd){
		this(beginInd, endInd, "", 1);
	}
	
	public String getPep(){
		//if(this.protein.length() > 0){
			//System.out.println("protein is: " + this.protein);
		//	return this.protein.substring((int)beginInd, (int)endInd);
			//return "";
		//}else{
			return this.fastaseq.getSubsequence(beginInd, endInd+1);
		//}
	}
	

	
	public String toString(){
		return this.getPep();
	}
}
