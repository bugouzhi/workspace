package Utils;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Perform logical operations on Sets
 * @author Jian
 *
 */
public class SetUtils {
	
	public static Set join(Set set1 , Set set2){
		Set joined = new HashSet(set1.size() + set2.size());
		joined.addAll(set1);
		joined.addAll(set2);
		return joined;
	}
	
	public static Set getIntersect(Set set1 , Set set2){
		Set intersect = new HashSet(set1.size() + set2.size());
		int count = 5;
		int i = 0;
		for(Iterator it = set1.iterator(); it.hasNext();){
			Object ele = it.next();
			if(set2.contains(ele)){
				intersect.add(ele);
			}
		}
		return intersect;
	}
	
	

}
