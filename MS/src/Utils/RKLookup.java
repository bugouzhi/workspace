package Utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.math.*;

/**
 * look up a set of strings using hashtable
 * @author Jian
 *
 */
public class RKLookup {
	private Map<Long, List<String>> table;
	private int prefixLen;
	int R = 256;
	long Q;
	long RLen;
	
	public RKLookup(Collection<String> patterns, int prefixLen){
		table = new HashMap<Long, List<String>>();
		this.Q = longRandomPrime();
		this.prefixLen = prefixLen;
		RLen = 1;
        for (int i = 1; i <= prefixLen-1; i++)
           RLen = (R * RLen) % Q;
		for(Iterator<String> it = patterns.iterator(); it.hasNext();){
			String pattern = it.next();
			//System.out.println("pattern " + pattern);
			if(pattern.length() < prefixLen){
				continue;
			}
			long key = hash(pattern, this.prefixLen);
			List<String> patternlist;
			if(this.table.containsKey(key)){
				patternlist = this.table.get(key);
			}else{
				patternlist = new ArrayList<String>();
			}
			patternlist.add(pattern);
			this.table.put(key, patternlist);
		}
		//System.out.println("Done creating table, keys: " + this.table.keySet().size() + "\tlooking for: " + patterns.size());
	}
	
	
	public Map<String, List<Integer>> matches(String toBeMatch){
		StringBuffer buffer = new StringBuffer(toBeMatch);
		Map<String, List<Integer>> matches = new HashMap<String, List<Integer>>();
		boolean match = false;
        int N = toBeMatch.length(); 
        if (N < this.prefixLen) return matches;
	    long txtHash = hash(toBeMatch, this.prefixLen); 
       // check for match at offset 0                         //still need to handle this case
        //if ((patHash == txtHash) && check(txt, 0))
        //    return 0;
       // check for hash match; if hash match, check for exact match
        
	    for (int i = this.prefixLen; i < N; i++) {
	            // Remove leading digit, add trailing digit, check for match. 
	            txtHash = (txtHash + Q - RLen*toBeMatch.charAt(i-prefixLen) % Q) % Q; 
	            txtHash = (txtHash*R + toBeMatch.charAt(i)) % Q; 
	            if(this.table.containsKey(txtHash)){
	            	List<String> patternlist = this.table.get(txtHash);
	            	for(int j = 0; j < patternlist.size(); j++){
	            		String pattern = patternlist.get(j);
	            		if(checkString(pattern, toBeMatch, i-prefixLen+1)){
	            			if(!matches.containsKey(pattern)){
	            				matches.put(pattern, new ArrayList());
	            			}
	            			List<Integer> pos = matches.get(pattern);
	            			pos.add(i-prefixLen+1);
	            		}
	            	}
	            }
	            //if(i % 1000000 == 0){
	            	//System.out.println(i);
	            //}
	    }
	    System.out.println("Done matching: " +"found matches:\t" + matches.keySet().size());
		  // no match
		  //  return N;
	    
	    return matches;
	}
	
	private boolean checkString(String pattern, String txt, int offset){
		for(int i = 0; i < pattern.length(); i++){
			if(i+offset >= txt.length()
					|| pattern.charAt(i) != txt.charAt(i+offset)){
				return false;
			}
		}
		return true;
	}
	
    // Compute hash for key[0..M-1]. 
    private long hash(String key, int M) { 
    	
        long h = 0; 
        for (int j = 0; j < M; j++) 
            h = (R * h + key.charAt(j)) % Q; 
        return h; 
    } 
    
    // a random 31-bit prime
    private static long longRandomPrime() {
        BigInteger prime = BigInteger.probablePrime(31, new Random());
        return prime.longValue();
    }

}
