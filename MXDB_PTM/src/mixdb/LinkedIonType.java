package mixdb;

public class LinkedIonType extends MixIonType{
	
	private boolean isLinked = false;
	
	public LinkedIonType(MixIonType type, boolean isLinked) {
		super(type);
		this.isLinked = isLinked;
	}

	public boolean isLinked() {
		return isLinked;
	}

	public void setLinked(boolean isLinked) {
		this.isLinked = isLinked;
	}
	
	public String toString(){
		if(isLinked)
			return super.toString()+"@linked";
		else
			return super.toString()+"@unlinked";
	}
	
}
