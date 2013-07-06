package mixdb;
/**
 * A simple peptide type 
 * @author Jian Wang
 *
 */
public class SimplePeptideType implements PeptideType{
	public static int SHORT = 0;
	public static int LONG = 1;
	private int pepCharge;
	private int length;
	
	public SimplePeptideType(int pepCharge, int pepLength){
		this.pepCharge = pepCharge;
		this.length = pepLength;
	}
	public int getPepCharge() {
		return pepCharge;
	}
	public void setPepCharge(int pepCharge) {
		this.pepCharge = pepCharge;
	}
	public int getLength() {
		return length;
	}
	public void setLength(int length) {
		this.length = length;
	}
	public String toString(){
		String len;
		if(this.length == SHORT){
			len = "short";
		}else{
			len = "long";
		}
		return this.pepCharge+"@"+len;
	}
	
	
}
