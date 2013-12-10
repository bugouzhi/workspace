package Utils;

public class StringUtils {
	public static String getStrippedSeq(String seq){
		return seq.replaceAll("[0-9\\+\\-\\.\\_]", "");
	}
	
	public static String getPepSeq(String seq){
		if(seq.charAt(1) == '.' && seq.charAt(seq.length()-2) == '.'){
			return seq.substring(2, seq.length()-2);
		}else{
			return seq;
		}
	}
}
