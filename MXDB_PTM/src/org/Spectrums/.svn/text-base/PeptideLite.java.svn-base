package org.Spectrums;
/**
 * A simple memory efficient implementation of peptides
 * @author Jian Wang
 *
 */
public class PeptideLite{
	private int beginInd;
	private int endInd;
	private int charge;
	private String protein;
	
	public PeptideLite(int beginInd, int endInd){
		this.beginInd = beginInd;
		this.endInd = endInd;
	}
	
	public PeptideLite(int beginInd, int endInd, String proteins, int charge){
		this.beginInd= beginInd;
		this.endInd = endInd;
		this.protein = proteins;
		this.charge = charge;
	}
	
	public String getPep(){
		return this.protein.substring(beginInd, endInd);
	}
}
