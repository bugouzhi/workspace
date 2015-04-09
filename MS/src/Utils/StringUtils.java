package Utils;

public class StringUtils {
	public static String getStrippedSeq(String seq){
		//System.out.println("before " + seq + "\t" + seq.replaceAll("[^A-Z]", ""));
		return seq.replaceAll("[^A-Z]", "");
	}
	
	
	public static String stripNeighborRes(String seq){
		if(seq.charAt(1) == '.' && seq.charAt(seq.length()-2) == '.'){
			return seq.substring(2, seq.length()-2);
		}else{
			return seq;
		}
	}
	/**
	 * Given a mixture result file from mixdb, reverse the stat values for
	 * the first and second peptide
	 * @param result
	 * @return
	 */
	public static String switchMixDBStat(String result){
		int[] indicesFrom = new int[]{3, 5, 7, 10, 12, 15, 16, 19, 20, 23, 25}; //indices to switch from;
		int[] indicesTo = new   int[]{4, 6, 8, 11, 13, 17, 18, 21, 22, 24, 26};
		String[] tokens = result.split("\\t");
		StringBuffer buff = new StringBuffer();
		for(int i = 0; i < indicesFrom.length; i ++){
			String tmp = tokens[indicesFrom[i]];
			tokens[indicesFrom[i]] = tokens[indicesTo[i]];
			tokens[indicesTo[i]] = tmp;
		}
		for(int i = 0; i < tokens.length; i++){
			buff.append(tokens[i]);
			if(i != tokens.length-1){
				buff.append("\t");
			}
		}
		return buff.toString();
	}
	
	public static int[] toIntArray(String[] arry){
		int[] values = new int[arry.length];
		for(int i = 0; i < arry.length; i++){
			values[i]  = Integer.parseInt(arry[i]);
		}
		return values;
	}
	
	public static double[] toDoubleArray(String[] arry){
		double[] values = new double[arry.length];
		for(int i = 0; i < arry.length; i++){
			values[i]  = Double.parseDouble(arry[i]);
		}
		return values;
	}
	
	
}
