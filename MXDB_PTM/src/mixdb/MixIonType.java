package mixdb;
/**
 * Represent mix ion type
 * @author Jian Wang
 *
 */
public class MixIonType extends SimpleIonType{

	public MixIonType(SimpleIonType type) {
		super(type);
	}
	
	
	public MixIonType(MixIonType type) {
		super(type);
		this.setPType(type.getPType());
	}
	
	public MixIonType(SimpleIonType type, int mixtureCharge, int Abundance) {
		super(type);
		SimplePeptideType pType = (SimplePeptideType) this.getPType();
		MixturePeptideType mType = new MixturePeptideType(pType.getPepCharge()
				, pType.getLength(), mixtureCharge, Abundance);
		this.setPType(mType);
		mType.setPeptideAbundance(Abundance);
	}
	
	public MixturePeptideType getPeptideType(){
		return (MixturePeptideType)this.getPType();
	}
		
	public int getAbundance(){
		return this.getPeptideType().getPeptideAbundance();
	}
	
	public boolean isHiAbundance(){
		return this.getPeptideType().getPeptideAbundance() == MixturePeptideType.HIGH_ABUNDANCE;
	}
	
	public int getMixtureCharge(){
		return ((MixturePeptideType)this.getPType()).getMixtureCharge();
	}
	
}
