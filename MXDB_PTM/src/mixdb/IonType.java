package mixdb;
/**
 * An interface represent iontype
 * @author Jian Wang
 *
 */
public interface IonType {
	public PeptideType getPType();
	//should consider migrate to this system, a should use a method to
	
	//compute the mass of an ions from the sequence masses
	//public double getIonMass(double prefixmass, double parentmass);
	
	//return the specific integer representation of this ion type
	//not sure should use it as a mehod cause hard to keep it consistent
	//or use a mapper to assign the index
	//public int getIonIndex()
}
