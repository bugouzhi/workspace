package mixdb;
/**
 * Represent co-eluting peptides in mixture spectra, distinguish
 * high-abundance peptides and low abundance peptide
 * @author Jian Wang
 *
 */
public class MixturePeptideType extends SimplePeptideType{
	public static int HIGH_ABUNDANCE = 0;
	public static int LOW_ABUNDANCE = 1;
	private int peptideAbundance;
	private int mixtureCharge;
	
	public MixturePeptideType(int pepCharge, int pepLength, int mixtureCharge, int abundance){
		super(pepCharge, pepLength);
		this.peptideAbundance=abundance;
		this.mixtureCharge = mixtureCharge;
	}
	
	public int getPeptideAbundance() {
		return peptideAbundance;
	}
	public void setPeptideAbundance(int peptideAbundance) {
		this.peptideAbundance = peptideAbundance;
	}
	
	
	public int getMixtureCharge() {
		return mixtureCharge;
	}

	public void setMixtureCharge(int mixtureCharge) {
		this.mixtureCharge = mixtureCharge;
	}

	public String toString(){
		if(this.peptideAbundance == HIGH_ABUNDANCE){
			return super.toString() + "@" + this.mixtureCharge + "@high";
		}else{
			return super.toString() + "@" + this.mixtureCharge + "@low";
		}
	}

}
