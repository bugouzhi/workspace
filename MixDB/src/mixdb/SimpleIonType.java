package mixdb;
/**
 * Represent the ion type generated from a single peptides
 * @author Jian Wang
 *
 */
public class SimpleIonType implements IonType{
		public static int PREFIX = 0;
		public static int SUFFIX = 1;
		private String type; //string representation of this ion type
		private double offset; //mass offset of this ion type
		private int charge;  //charge of this ion type
		private int direction;
		private PeptideType pType;
		
		
		public SimpleIonType(String type, double offset, int charge, int direction, int pepCharge, int pepLength){
			this(type, offset, charge, direction, new SimplePeptideType(pepCharge, pepLength));
		}
		
		public SimpleIonType(String type, double offset, int charge, int direction, PeptideType pType){
			this.type = type;
			this.offset = offset;
			this.charge = charge;
			this.direction = direction;
			this.pType = pType;
		}
		
		/**
		 * copy constructor
		 * @param type
		 */
		public SimpleIonType(SimpleIonType type){
			this(type.getType(), type.getOffset(), type.getCharge(), 
					type.getDirection(), type.getPepCharge(), 
					((SimplePeptideType)type.getPType()).getLength());
		}
		
		public int getPepCharge() {
			return ((SimplePeptideType)pType).getPepCharge();
		}
		
		public int getPepLength() {
			return ((SimplePeptideType)pType).getLength();
		}
		
		public String getType() {
			return type;
		}
		public void setType(String type) {
			this.type = type;
		}
		public double getOffset() {
			return offset;
		}
		public PeptideType getPType() {
			return pType;
		}

		public void setPType(PeptideType type) {
			pType = type;
		}

		public void setOffset(double offset) {
			this.offset = offset;
		}
		public int getCharge() {
			return charge;
		}
		public void setCharge(int charge) {
			this.charge = charge;
		}
		
		public String toString(){
			return this.type+"@"+this.charge+"@"+this.pType.toString();
		}
		
		public int getDirection() {
			return direction;
		}

		public void setDirection(int direction) {
			this.direction = direction;
		}
		
		
}
