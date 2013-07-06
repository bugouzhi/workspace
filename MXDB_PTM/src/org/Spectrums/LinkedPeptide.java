package org.Spectrums;

import java.util.List;

public class LinkedPeptide extends Peptide{
	//private String crossLinker = null;
	//private double crossLinkerMass = 0;
	Peptide[] peptides; //treated linked peptide as two peptide with PTM

	public LinkedPeptide(){
		
	}
	
	public LinkedPeptide(String pep, int charge){
		String[] peps = pep.split("--");
		int pos1 = peps[0].indexOf('K')+1;
		int pos2 = peps[1].indexOf('K')+1;
		createLinkedPeptide(pep, charge, pos1, pos2);
	}
	
	public LinkedPeptide(String pep, int charge, int position1, int position2){
		createLinkedPeptide(pep, charge, position1, position2);
	}
	
	public Peptide[] getPeptides(){
		return this.peptides;
	}
	private void createLinkedPeptide(String pep, int charge, int position1, int position2){
		String[] peps = pep.split("--");
		Peptide p1 = new Peptide(peps[0], 1);
		Peptide p2 = new Peptide(peps[1], 1);
		
		this.peptides = new Peptide[2];
		this.peptides[0] = p1;
		this.peptides[1] = p2;
		double mass = (p1.getParentmass() + p2.getParentmass()
				+ Mass.DSSLINKER_MASS 
				+ Mass.PROTON_MASS*(charge-2))/charge;	
		double massShift1 = p1.getParentmass();
		double massShift2 = p2.getParentmass();
		p1.insertPTM(position1, massShift2+Mass.DSSLINKER_MASS-Mass.PROTON_MASS);
		p2.insertPTM(position2, massShift1+Mass.DSSLINKER_MASS-Mass.PROTON_MASS);
		p1.setLinkedPos(position1);
		p2.setLinkedPos(position2);
//		System.out.println("we have ptms: " + p1.getPos().length);
//		System.out.println("we have ptms: " + p2.getPos().length);
		this.setCharge((short)charge);
//		System.out.println("peptides is " + pep);
		this.setPeptide(p1 + "--" + p2);
//		double mass = (PeptideMassAnalysis.computeMolecularMass(peps[0])
//				+ PeptideMassAnalysis.computeMolecularMass(peps[1]) 
//				+ Mass.DSSLINKER_MASS 
//				+ Mass.PROTON_MASS*(charge))/charge;	
		this.setParentmass(mass);
	}
	
	public LinkedPeptide(Peptide pep1, Peptide pep2, int charge, int position1, int position2){
		Peptide p1 = new Peptide(pep1);
		Peptide p2 = new Peptide(pep2);
		this.peptides = new Peptide[2];
		this.peptides[0] = p1;
		this.peptides[1] = p2;
		double mass = (p1.getParentmass() + p2.getParentmass()
				+ Mass.DSSLINKER_MASS 
				+ Mass.PROTON_MASS*(charge-2))/charge;	
		double massShift1 = p1.getParentmass();
		double massShift2 = p2.getParentmass();
 		p1.insertPTM(position1, massShift2+Mass.DSSLINKER_MASS-Mass.PROTON_MASS);
		p2.insertPTM(position2, massShift1+Mass.DSSLINKER_MASS-Mass.PROTON_MASS);
		p1.setLinkedPos(position1);
		p2.setLinkedPos(position2);
		//System.out.println("we have ptms: " + p1.getPos().length);
		//System.out.println("we have ptms: " + p2.getPos().length);
		this.setCharge((short)charge);
		//this.setPeptide(p1.getPeptide() + "--" + p2.getPeptide());	
		this.setParentmass(mass);
	}
	
	
	public LinkedPeptide(Peptide p1, Peptide p2, int charge){
		this.peptides = new Peptide[2];
		this.peptides[0] = p1;
		this.peptides[1] = p2;
		
		double mass1 = p1.getParentmass()*p1.getCharge() - Mass.PROTON_MASS*p1.getCharge();
		double mass2 = p2.getParentmass()*p2.getCharge() - Mass.PROTON_MASS*p2.getCharge();	
		
		double pMass1 = -1;
		double pMass2 = -2;
		int index1 = -1;
		int index2 = -1;
		//recompute theoretical masses, since it may not have to 
		//correct from the half-linked peptide
		//System.out.println(p1.getPos()[0] +"\t" + p1.getLinkedPos() + "\t" + p2.getPos()[0] + "\t" + p2.getLinkedPos());
		double[] ptmMasses1 = p1.getPtmmasses();
		double[] ptmMasses2 = p2.getPtmmasses();
		for(int i = 0; i < ptmMasses1.length; i++){
			if(p1.getPos()[i]== p1.getLinkedPos()){
				pMass1= mass1 - p1.getPtmmasses()[i];
				index1 = i;
			}
		}
		
		for(int i = 0; i < ptmMasses2.length; i++){
			if(p2.getPos()[i] == p2.getLinkedPos()){
				pMass2 = mass2 - p2.getPtmmasses()[i];
				index2 = i;
			}
		}
		ptmMasses1[index1]=pMass2+Mass.DSSLINKER_MASS;
		ptmMasses2[index2]=pMass1+Mass.DSSLINKER_MASS;
		p1.setPtmmasses(ptmMasses1);
		p2.setPtmmasses(ptmMasses2);
		double mass = (pMass1 + pMass2
				+ Mass.DSSLINKER_MASS 
				+ Mass.PROTON_MASS*(charge))/charge;
//		//System.out.println("we have ptms: " + p1.getPos().length);
//		//System.out.println("we have ptms: " + p2.getPos().length);
		this.setCharge((short)charge);
		this.setPeptide(p1.getPeptide() + "--" + p2.getPeptide());	
		this.setParentmass(mass);
	}

	public static LinkedPeptide createLinkedPeptide(String pep, int charge){
		String[] peps = pep.split("--");
		int position1 = peps[0].indexOf('+');
		int position2 = peps[1].indexOf('+');
		peps[0] = peps[0].replaceAll("[0-9\\.\\+]", "");
		peps[1] = peps[1].replaceAll("[0-9\\.\\+]", "");
		return new LinkedPeptide(peps[0] + "--" + peps[1], charge, position1, position2);
	}
	
	public LinkedPeptide(LinkedPeptide lp){
		Peptide p1 = new Peptide(lp.peptides[0]);
		Peptide p2 = new Peptide(lp.peptides[1]);
		this.peptides = new Peptide[]{p1, p2};
		this.setParentmass(lp.getParentmass());
		this.setCharge(lp.getCharge());
		//this.setPeptide(p1.getPeptide() + "--" + p2.getPeptide());
		this.setPeptide(lp.getPeptide());
	}
	
	public static int getMinLinkedCharge(int linkedCharge){
		if(linkedCharge == 1){
			return 1;
		}else if(linkedCharge == 2){
			return 1;
		}else if(linkedCharge == 3){
			return 2;
		}else if(linkedCharge == 4){
			return 2;
		}else if(linkedCharge == 5){
			return 3;
		}else{
			return 4;
		}
			
	}
	
	public static int getMaxLinkedCharge(int linkedCharge){
		if(linkedCharge == 1){
			return 1;
		}else if(linkedCharge == 2){
			return 2;
		}else if(linkedCharge == 3){
			return 2;
		}else if(linkedCharge == 4){
			return 3;
		}else if(linkedCharge == 5){
			return 4;
		}else{
			return 5;
		}
	}
	
	public String toString(){
		return peptides[0] + "--" + peptides[1];
	}
	
	public static int transformPeakCharge(int peakCharge, int pepCharge){
		if(pepCharge == 2){
			return peakCharge;
		}else if(pepCharge == 3){
			return peakCharge - 1;
		}else if(pepCharge == 4){
			return peakCharge - 1;
		}else if(pepCharge == 5){
			return peakCharge - 2;
		}else{
			if(pepCharge > 5){
				return peakCharge-3;
			}else{
				return peakCharge - 2;
			}
		}
	}
	public static void testLinkedPeptideMass(){
		String filename = "..\\mixture_linked\\t";
		List<String> lines = Utils.FileIOUtils.createListFromFile(filename);
		for(int i = 0; i < lines.size(); i++){
			String[] tokens = lines.get(i).split("\\s+");
			double mass = Double.parseDouble(tokens[3]);
			int charge = 3;
			if(mass < 550){
				charge =4;
			}
			LinkedPeptide lp = new LinkedPeptide(tokens[4] + "--" + tokens[6],charge);
			double pmass = lp.getParentmass();
			double pmass12 = pmass +(-138.0680 + 150.1439)/charge;
			//double diff = Math.abs(pmass - mass + pmass/100000);
			//double diff2 = Math.abs(pmass12 - mass + pmass/10000);
			double diff = Math.abs((pmass - mass)*1000000/pmass+100);
			double diff2 = Math.abs((pmass12 - mass)*1000000/pmass12+100);
			System.out.println(tokens[0] + " " + tokens[1] + " " + tokens[2] + "\t" + 
					lp + "\t" +  pmass + "\t" + pmass12 + "\t" +  mass + "\t" 
					+ diff + "\t" + diff2);
		}
	}
	
	public static void main(String[] args){
		testLinkedPeptideMass();
	}
}
