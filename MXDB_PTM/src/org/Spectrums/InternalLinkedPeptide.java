package org.Spectrums;

public class InternalLinkedPeptide extends LinkedPeptide{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1324234324L;
	
	
	public InternalLinkedPeptide(String pep, int charge){
		super(pep, charge);
	}
	
	public InternalLinkedPeptide(Peptide p, int charge, int position1, int position2){
		//System.out.println("creating internal linked: " + p.getPeptide() + "\t" + charge + "\t" + position1 + "\t" + position2);
		String pep = p.getPeptide();
		this.setPeptide(pep);
		Peptide p1 = new Peptide(pep.substring(0, position1), 1);
		//System.out.println("peptide1: " + p1.getPeptide());
		Peptide p2 = new Peptide(pep.substring(position2-1), 1);
		//System.out.println("peptide2: " + p2.getPeptide());
		this.peptides = new Peptide[2];
		this.peptides[0] = p1;
		this.peptides[1] = p2;
		double mass = ((p.getParentmass() - Mass.PROTON_MASS*p.getCharge())
				+ Mass.DSSLINKER_MASS); 
		double massShift1 = mass - (p1.getParentmass() - p1.getCharge()*Mass.PROTON_MASS - Mass.WATER);
		double massShift2 = mass - (p2.getParentmass() - p2.getCharge()*Mass.PROTON_MASS - Mass.WATER);
 		p1.insertPTM(position1, massShift1);
		p2.insertPTM(1, massShift2);
		p1.setLinkedPos(position1);
		p2.setLinkedPos(1);
//		System.out.println("we have ptms: " + p1.getPos().length);
//		System.out.println("we have ptms: " + p2.getPos().length);
		this.setCharge((short)charge);
//		System.out.println("peptides is " + pep);
		this.setPeptide(p1  + "--[" + pep.substring(position1, position2-1) +"]-" + p2);
		mass = (mass + Mass.PROTON_MASS*charge)/charge;
		this.setParentmass(mass);
	}
	
	public String toString(){
		return this.getPeptide();
	}
	
	
}
