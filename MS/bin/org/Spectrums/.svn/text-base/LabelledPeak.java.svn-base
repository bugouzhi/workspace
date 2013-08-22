package org.Spectrums;
/**
 * A peak with annotation
 * @author jian wang
 *
 */
public class LabelledPeak extends Peak{
	private String type;
	private short pos;
	private short charge;
	private Peptide pep;
	public static double DEFAULT_INTENS = 1000;
		
	public LabelledPeak(double moz, double intensity, String type, short pos, short charge){
		super(moz, intensity);
		this.type = type;
		this.pos = pos;
		this.charge = charge;
	}
	
	public LabelledPeak transferLabel(Peak p){
		return new LabelledPeak(p.getMass(), p.getIntensity(), this.getType(), this.getPos(), this.getCharge());
	}
	
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public short getPos() {
		return pos;
	}
	public void setPos(short pos) {
		this.pos = pos;
	}
	public short getCharge() {
		return charge;
	}
	
	public void setCharge(short charge) {
		this.charge = charge;
	}
	
	public String toString(){
		if(this.pep == null){
			return "\t" + this.printType()+ "\t" + this.pos+"@"+this.charge+"\t" + this.getMass() + "\t"+this.getIntensity();

		}
		return this.pep.getPeptide() + ":\t" + this.printType()+ "\t" + this.pos+"@"+this.charge+"\t" + this.getMass() + "\t"+this.getIntensity();
	}
	
	public String printType(){
		String formatType = this.type.charAt(0) + "\t" + this.type;
		return formatType;
	}
	public Peptide getPep() {
		return pep;
	}

	public void setPep(Peptide p) {
		this.pep = p;
	}
	
	public boolean isPrefixPeak(){
		return this.type.contains("b") || this.type.contains("a");
	}
	
	public boolean isSuffixPeak(){
		return this.type.contains("y");
	}
	
	public boolean isPrimary(){
		return this.type.equals("b") || this.type.equals("y") || this.type.equals("a") || this.type.equals("b(X)") || this.type.equals("y(X)");
			//|| this.type.equals("b(iso)") || this.type.equals("y(iso)");
	}
	public String getPeakType(){
		return this.getType()+"@"+this.getPos()+"@"+this.getCharge();
	}

}
