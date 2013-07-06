package org.Spectrums;

public class MixturePeak extends LabelledPeak{

	/**
	 * Serial ID for printing object
	 */
	private static final long serialVersionUID = 5595832889800794039L;
	private int peptideIndex;
	private Spectrum parent;
	public Spectrum getParent() {
		return parent;
	}

	public void setParent(Spectrum parent) {
		this.parent = parent;
	}

	public MixturePeak(double moz, double intensity, String type, short pos, short charge, int peptideIndex){
		super(moz, intensity, type, pos, charge);
		this.peptideIndex =  peptideIndex;
	}
	
	public MixturePeak(LabelledPeak p, int peptideIndex){
		this(p.getMass(), p.getIntensity(), p.getType(), p.getPos(), p.getCharge(), peptideIndex);
		this.setPep(p.getPep());
	}
	
	public int getPeptideIndex() {
		return peptideIndex;
	}
	
	public void setPeptideIndex(int peptideIndex) {
		this.peptideIndex = peptideIndex;
	}
	
	public String toString(){
		if(this.getPep() == null){
			return "\t" + this.printType()+ "\t" + this.getPos()+"@"+this.getCharge()
			+"@"+this.getPep().getCharge()+"@"+ this.getParent().charge
			+"\t" + this.getMass() + "\t"+this.getIntensity() + "\t" + getPeptideIndex();	

		}
		return this.getPep().getPeptide() + ":\t" + this.printType()+ "\t" + this.getPos()+"@"+this.getCharge()
			+"@"+this.getPep().getCharge()+"@"+ this.getParent().charge
			+"\t" + this.getMass() + "\t"+this.getIntensity() + "\t" + getPeptideIndex();
	}
	
}
