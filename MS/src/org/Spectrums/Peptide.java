package org.Spectrums;
import java.io.Serializable;
import java.util.List;
import java.util.ArrayList;

import org.apache.commons.collections.MultiMap;

/**
 * Represent a peptide 
 * @author jian wang
 *
 */
public class Peptide implements Serializable{
	public static final long serialVersionUID = 12387183721L;
	private static int[] EMPTY_POS = new int[0];
	private static double[] EMPTY_PTMS = new double[0];
	private String peptide;
	private char[] pepseq;
	private int beginIndex;
	private int endIndex;
	private int[] pos = EMPTY_POS;
	private double[] ptmmasses = EMPTY_PTMS;
	private short charge;
	private double parentmass;
	private boolean isDecoy = false;
	private int linkedPos = -1; //-1 can mean it has not been setted yet
	
	public int getLinkedPos() {
		return linkedPos;
	}

	public void setLinkedPos(int linkedPos) {
		this.linkedPos = linkedPos;
	}

	public boolean isDecoy() {
		return isDecoy;
	}

	public void setDecoy(boolean isDecoy) {
		this.isDecoy = isDecoy;
	}

	public double getParentmass() {
		if(this.parentmass == 0){
			this.parentmass = PeptideMassAnalysis.computeMbyZ(this, this.charge);
		}
		return this.parentmass;
	}
	
	public void setParentmass(double parentmass) {
		this.parentmass = parentmass;
	}
	
	public Peptide(){
		
	}
	
	public Peptide(String peptide){
		if(peptide.startsWith("r") ){
			peptide = peptide.substring(1);
		}
		if(peptide.startsWith("X_") ){
			peptide = peptide.substring(2);
		}

		this.peptide = parsePeptide2(peptide);
	}
	
//	public Peptide(String pep, int charge){
//		this.peptide = pep;
//		this.charge = (short)charge;
//	}
	
	public Peptide(String pep, int charge){
		//System.out.println("creating peptide: " + pep);
		this.charge = (short)charge;
		this.peptide = parsePeptide2(pep+"."+charge);
		//this.peptide = pep;
	}
	
	//copy constructor
	public Peptide(Peptide p){
		this();
		this.peptide = p.getPeptide();
		this.charge = p.getCharge();
		this.parentmass = p.parentmass;
		this.pos = new int[p.pos.length];
		this.ptmmasses =  new double[p.ptmmasses.length];
		this.linkedPos = p.linkedPos;
		for(int i = 0; i < p.pos.length; i++){
			this.pos[i] = p.pos[i];
		}
		for(int i = 0; i < p.ptmmasses.length; i++){
			this.ptmmasses[i] = p.ptmmasses[i];
		}
	}
	
	//@TODO needed to update masses when modify the peptide
	public Peptide reverse(){
		Peptide p = new Peptide(this);
		p.peptide = new StringBuffer(this.peptide).reverse().toString();
		return p;
	}
	
	public static boolean isValidPeptide(String pep){
		return !pep.contains("Scan");             //cheap way of checking peptide, will need to extends later
	}
	
	private String parsePeptide(String pep){
		//System.out.println("peptide is: " + pep);
		if(pep.contains("r")){
			parsePeptide(pep.substring(1));
		}
		if(pep.contains("&")){
			parseMultiplePeptide(pep);
		}
		String[] tokens = pep.split("\\.");
		int chargeInd = pep.lastIndexOf('.');
		this.charge = Short.parseShort(pep.substring(chargeInd+1));
		pep = pep.substring(0, chargeInd);
		System.out.println("peptide is: " + pep);
		StringBuffer raw = new StringBuffer(), ptm = new StringBuffer();
		for(int i = 0; i < pep.length(); i++){
			if(pep.charAt(i) == '+' || pep.charAt(i) == '-' || pep.charAt(i) == '.'
				|| Character.isDigit(pep.charAt(i))){
				ptm.append(pep.charAt(i));
			}else{
				raw.append(pep.charAt(i));
			}
		}
		if(ptm.length() > 2){
			ptmmasses = new double[]{Double.parseDouble(ptm.toString())};
			pos = new int[]{pep.indexOf('+')-1};
			System.out.println("ptm: " + " +"+ ptmmasses[0] + "@"+pos[0]);
			//getAccuratePTMmass();
		}
		return raw.toString();
	}
	
	private String parsePeptide2(String pep){
		//System.out.println("peptide is: " + pep);
		if(pep.contains("r")){
			parsePeptide(pep.substring(1));
		}
		if(pep.contains("&")){
			parseMultiplePeptide(pep);
		}
		String[] tokens = pep.split("\\.");
		int chargeInd = pep.lastIndexOf('.');
		this.charge = Short.parseShort(pep.substring(chargeInd+1));
		pep = pep.substring(0, chargeInd);
		StringBuffer raw = new StringBuffer(), ptm = new StringBuffer();
		List<Double> masses = new ArrayList<Double>();
		List<Integer> positions = new ArrayList<Integer>();
		int resPos = 0, ptmPos = 0;
		boolean isSpectraSTFormat=false;
		for(int i = 0; i < pep.length(); i++){
			if(pep.charAt(i) == 'n'){
				continue;
			}
			if(pep.charAt(i) == '+' || pep.charAt(i) == '-' || pep.charAt(i) == '['){
				if(ptm.length() > 0){					                                                                                
					if(isSpectraSTFormat){
						masses.add(Double.parseDouble(ptm.toString())- Mass.getAAMass(raw.charAt(ptmPos-1))); //in spectraST format, inside the [] bracket
						isSpectraSTFormat = false;                                                            //it has residue+ptm rather than just the ptm offset
					}else{
						masses.add(Double.parseDouble(ptm.toString()));
					}
					positions.add(new Integer(ptmPos));
					ptm = new StringBuffer();
					
				}
				if(pep.charAt(i) == '-'){
					ptm.append('-');
				}
				if(pep.charAt(i+1) == '-'){
					ptm.append('-');
					i++;
				}
				if(pep.charAt(i) == '['){
					isSpectraSTFormat = true;
				}
				
				ptmPos = resPos;
				if(ptmPos == 0){   //not differetiating n-term of mod on first residue
					ptmPos = 1;    //makes all mod on n-term as if they are on first residue
				}
			}else if(pep.charAt(i) == ']'){
				
			}else if(Character.isDigit(pep.charAt(i)) || pep.charAt(i) == '.'){
				ptm.append(pep.charAt(i));
			}else{
				raw.append(pep.charAt(i));
				resPos++;
			}
		}
		if(ptm.length() > 0){ //taking care of last ptm
			if(isSpectraSTFormat){
				masses.add(Double.parseDouble(ptm.toString())- Mass.getAAMass(raw.charAt(ptmPos-1))); //in spectraST format, inside the [] bracket
				isSpectraSTFormat = false;                                                            //it has residue+ptm rather than just the ptm offset
			}else{
				masses.add(Double.parseDouble(ptm.toString()));
			}
			positions.add(new Integer(ptmPos));
			ptm = new StringBuffer();
		}
		if(positions.size() > 0){
			this.ptmmasses = new double[masses.size()];
			this.pos = new int[positions.size()];
			for(int i = 0; i < masses.size(); i++){
				this.ptmmasses[i] = masses.get(i);
				this.pos[i] = positions.get(i);
				//System.out.println("ptm: " +  ptmmasses[i] + "@"+pos[i]);
				//getAccuratePTMmass();
			}
		}
		return raw.toString();
	}
	
	
	private String parseMultiplePeptide(String pep){
		return "";
	}
	
	/**\
	 * Obtain more accurate ptmmasses from table
	 * @param ptmmasses
	 */
	//for now just for DSP linkers
	private void getAccuratePTMmass(){
		if(this.ptmmasses[0] == 145){
			this.ptmmasses[0] = Mass.DSPDANGLE_MASS;
		}
	}
	
	public boolean hasPTMs(){
		return this.pos != null;
	}

	public String getPeptide() {
		return this.peptide;
	}

	public void setPeptide(String peptide) {
		this.peptide = peptide;
		if(!peptide.contains("--")){
			this.parentmass = PeptideMassAnalysis.computeMbyZ(this, this.getCharge());			
		}
	}

	public int[] getPos() {
		return pos;
	}

	public void setPos(int[] pos) {
		this.pos = pos;
	}

	public double[] getPtmmasses() {
		return ptmmasses;
	}

	public void setPtmmasses(double[] ptmmasses) {
		this.ptmmasses = ptmmasses;
		this.parentmass = PeptideMassAnalysis.computeMbyZ(this, this.getCharge());
	}

	public short getCharge() {
		return charge;
	}

	public void setCharge(short charge) {
		double uncharged = this.parentmass * this.charge - this.charge*Mass.PROTON_MASS;
		this.charge = charge;
		this.parentmass = (uncharged + this.charge*Mass.PROTON_MASS)/this.charge;
	}
	//does index for ptm needed to be in order?
	public void insertPTM(int pos, double mass){
		//System.out.println("before we have: " + this.pos.length + "\t" + this.ptmmasses.length);
		if(this.pos == null){
			this.pos = new int[]{pos};
			this.ptmmasses = new double[]{mass};
		}else{
			int[] newpos = new int[this.pos.length+1];
			double[] newmasses = new double[this.ptmmasses.length+1];
			int j = 0;
			boolean hasInsert = false;
			for(int i = 0; i < newpos.length; i++){
				if(j < this.pos.length && this.pos[j] < pos){
					newpos[i]=this.pos[j];
					newmasses[i] = this.ptmmasses[j++];
				}else if(!hasInsert){
					newpos[i] = pos;
					newmasses[i] = mass;
					hasInsert = true;
				}else{
					newpos[i]=this.pos[j];
					newmasses[i]=this.ptmmasses[j++];
				}
			}
//			for(;j < this.pos.length; j++){
//				newpos[j] = this.pos[j];
//				newmasses[j]=this.ptmmasses[j];
//			}
//			newpos[j] = pos;
//			newmasses[j] = mass;
			this.pos = newpos;
			this.ptmmasses = newmasses;
		}
		this.parentmass =  PeptideMassAnalysis.computeMbyZ(this, this.charge);
		if(this.pos.length != this.ptmmasses.length){
			System.out.println("warining ptm is not created consistently");
		}
	}
	
	public static List<Peptide> insertPTM(List<Peptide> origList, double mass, char[] residues, int maxPTM){
		List<Peptide> modList = new ArrayList(origList.size());
		//System.out.println("we start with " + origList.size());
		for(int i = 0; i < origList.size(); i++){
			Peptide p = origList.get(i);
			for(int j=0; j < residues.length; j++){
				int index = p.getPeptide().indexOf(residues[j]);
				while(index >= 0){
					if(p.getPtmmasses().length < maxPTM){
						Peptide variant = new Peptide(p);
						//Peptide variant = p;
						variant.insertPTM(index+1, mass);
						modList.add(variant);
					}
					index = p.getPeptide().indexOf(residues[j], index+1);
				}
			}			
		}
		return modList;
	}
	
	public static List<Peptide> insertPTM(List<Peptide> origList, double mass, int position, int maxPTM){
		List<Peptide> modList = new ArrayList(origList.size());
		//System.out.println("we start with " + origList.size());
		for(int i = 0; i < origList.size(); i++){
			Peptide p = origList.get(i);
			if(p.getPeptide().length() >= position && p.getPtmmasses().length < maxPTM){
				Peptide variant = new Peptide(p);
				variant.insertPTM(position, mass);
				modList.add(variant);
			}
		}
		return modList;
	}
	
	public void createDSPLinkerPTM(){
		createDSPLinkerPTM(new int[] {this.peptide.indexOf('K')});
	}
	
	public void createDSPLinkerPTM(int[] pos){
		//System.out.println("setting DSP at position: " + pos[0]);
		this.pos = pos;
		this.ptmmasses = new double[]{145}; 
	}
	
	//check whether the peptide sequence is the same
	//this comparison ignore charge state as well as
	//ohter information such as ptm etc...
	public boolean isSamePeptideSeq(Peptide p){
		return this.peptide.equals(p.getPeptide());
	}
	
	public List<Integer> getLysPositions(){
		List<Integer> positions = new ArrayList();
		int i = this.peptide.indexOf('K');
		while(i >= 0){
			positions.add(new Integer(i));
			i = this.peptide.indexOf('K', i+1);
		}
		return positions;
	}
	
	public String toString(){
		StringBuffer buff = new StringBuffer(this.peptide);
		if(this.isDecoy){
			//buff.insert(0, 'r');
		}
		int offset=0;
		for(int i = 0; i < this.pos.length; i++){
			String ptm;
			if(ptmmasses[i]>0){
				ptm = String.format("+%.3f", this.ptmmasses[i]);
			}else{
				ptm = String.format("%.3f", this.ptmmasses[i]);
			}
			buff.insert(this.pos[i]+offset, ptm);
			offset+=ptm.length();
		}
		return buff.toString();
	}
}
