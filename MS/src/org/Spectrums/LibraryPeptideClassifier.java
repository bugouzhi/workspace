package org.Spectrums;
/**
 * Classify peptides according to definition from NIST 
 * The classes are used to help improve quality of spectral library
 * @author Jian
 *
 */
public class LibraryPeptideClassifier {

	
	
	
	public interface PeptideClass{
		public boolean isClass(String peptide);
	}
	
	public class PeptideClass1{
		public boolean isClass(String peptide){
			return true;
		}
	}
	
	public class PeptideClass2{
		public boolean isClass(String peptide){
			return true;
		}
	}
	public class PeptideClass3{
		public boolean isClass(String peptide){
			return true;
		}
	}
	public class PeptideClass4{
		public boolean isClass(String peptide){
			return true;
		}
	}
	public class PeptideClass5{
		public boolean isClass(String peptide){
			return true;
		}
	}
	public class PeptideClass6{
		public boolean isClass(String peptide){
			return true;
		}
	}
	public class PeptideClass7{
		public boolean isClass(String peptide){
			return true;
		}
	}
	public class PeptideClass8{
		public boolean isClass(String peptide){
			return true; 
		}
	}
	public class PeptideClass9{
		public boolean isClass(String peptide){
			return true;
		}
	}
	
	
}
