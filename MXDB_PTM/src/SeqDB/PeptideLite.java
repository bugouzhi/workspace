package SeqDB;
import java.io.Serializable;
import sequences.FastaSequence;

public class PeptideLite implements Serializable{

	/**
	 * A simple memory efficient implementation of peptides
	 * @author Jian Wang
	 *
	 */
		private static final long serialVersionUID = 230483204832L;
		private static String EMPTYSTR = "";
		private int beginInd;
		private int endInd;
		private int parentmass;
		//private int charge;
		//private String protein;
		//private FastaSequence fastaseq;
		
//		}
		
		public PeptideLite(int beginInd, int endInd, int parentMass){
			this.beginInd = beginInd;
			this.endInd = this.endInd;
			this.parentmass = parentMass;
		}
		
		public int getParentMass(){
			return this.parentmass;
		}

}
