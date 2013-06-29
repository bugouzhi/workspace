package mixdb;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
/**
 * Generate a mapping between ion type and an index
 * this index is used to represent ion type in a theoretical spectrum
 * @author Jian Wang
 *
 */
public class IonTypeMapper {
	private Map<Integer, IonType> mapper;
	private Map<IonType, Integer> mapper2;
	private int currentIndex=0;
	public IonTypeMapper(){
		this.mapper = new HashMap();
		this.mapper2 = new HashMap();

	}
	
	public void addIonType(IonType newType){
		Integer index = new Integer(currentIndex);
		this.mapper.put(index, newType);
		this.mapper2.put(newType, index);
		currentIndex++;
	}
	/**
	 * return the index for this ion type, return negative one 
	 * if iontype is not found
	 * @param type
	 * @return
	 */
	public int getIndex(IonType type){
		if(this.mapper2.containsKey(type)){
			return this.mapper2.get(type).intValue();
		}else{
			throw new IllegalArgumentException("ion type not valid: " + type);
		}
	}
	
	public IonType getIonType(int index){
		Integer ind = new Integer(index);
		if(this.mapper.containsKey(ind)){
			return this.mapper.get(ind);
		}else{
			throw new IllegalArgumentException("ion type index not valid: " + index);
		}
	}
	
	public int getSize(){
		return currentIndex;
	}
	
	public static IonTypeMapper createIonTypeMap(Collection<IonType> ionTypes){
		IonTypeMapper ionTypeMap = new IonTypeMapper();
		for(Iterator<IonType> it = ionTypes.iterator(); it.hasNext();){
			ionTypeMap.addIonType(it.next());
		}
		return ionTypeMap;
	}
	
		
}
