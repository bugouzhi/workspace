package org.Spectrums;



/**
 * operator for a individual element of lookup table
 * @author Jian Wang
 *
 */
public interface TableElementOperator {
	public double eval(double val1, double val2);
	public static class Divider implements TableElementOperator{
		public static Divider d = new Divider();
		private Divider(){
			
		}
		public double eval(double val1, double val2) {
			return val1/val2;
		}
		
	}
	
	public static class Log implements TableElementOperator{
		public static Log l = new Log();
		private Log(){
			
		}
		public double eval(double val1, double val2) {
			return Math.log(val1);
		}
		
	}
}
