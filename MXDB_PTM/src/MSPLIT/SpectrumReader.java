package MSPLIT;
/**
 * Read in spectrum from various format
 * @author jian wang
 *
 */
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Vector;


public class SpectrumReader {
	public static List<Spectrum> readSpectrumFromMGF(String file){
		List<Spectrum> v = new Vector<Spectrum>();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			Spectrum s = readSpectrumFromMGF(bf);
			while(s != null){
				v.add(s);
				s = readSpectrumFromMGF(bf);
			}		 
			
		}catch(IOException ioe){
			System.out.println("Cannot Open MGF file");
			System.out.println(ioe.getMessage());
		}
		return v;
	}
	
	public static List<Spectrum> readSpectrumFromMSP(String file){
		List<Spectrum> v = new Vector<Spectrum>();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			Spectrum s = readSpectrumFromMSP(bf);
			while(s != null){
				v.add(s);
				s = readSpectrumFromMSP(bf);
			}		 
			
		}catch(IOException ioe){
			System.out.println("Cannot Open MGF file");
			System.out.println(ioe.getMessage());
		}
		return v;
	}
	
	private static Spectrum readSpectrumFromMGF(BufferedReader bf) {
		SimpleSpectrum s = new SimpleSpectrum();
		try{
			String line;
			do{
				line = bf.readLine();
				//System.out.println(line);
			}while(line != null && !line.equals("BEGIN IONS"));
			// a flag tell us whether we are in peak sections
            //in the spectrum file
			boolean isPeaks = false; 
			String[] token;
			String charge;
			line = bf.readLine();
			while(line != null && !line.equals("END IONS")){
				//System.out.println("line is " + line);
				if(line.startsWith("PEPMASS")){
					s.setParentMass(Double.valueOf(((line.split("="))[1])));
					//isPeaks = true;  //for real data last line of mgf is pepmass
				}else if(line.startsWith("CHARGE")){
					charge = ((line.split("="))[1]);
					if(charge.startsWith("+")){
						charge = charge.substring(1);
					}
					s.setCharge(Integer.valueOf(charge));
					
				}else if(line.startsWith("PEPSEQ")){
					//this.peptide = (line.split("="))[1];
					//this.peptide = this.peptide + "." + this.charge;
				}else if(line.startsWith("TITLE")){
					//this.peptide = (line.split("="))[1];
					//this.peptide = this.peptide + "." + this.charge;
					s.setSpectrumName((line.split("="))[1]);
				}else if(line.startsWith("PEPMOD")){
					//just put some arbitary mod for now
					//mainly use this as a flag to tell spectra with mod
					//from those without any mod
					//this.modMass = 50;
					//this.modPos = 3;
					
				}else if(Character.isDigit(line.charAt(0))){
					//System.out.println("currelint is: " + line + "\tfirst: " + line.charAt(0));
					isPeaks = true;
					
				}
				if(isPeaks){
					token = line.split("\\s+");
					s.addPeak(new SimplePeak(Double.valueOf(token[0]) ,//* 0.9995, //note multiply by 0.9995 to make peaks center around integer values
						  Double.valueOf(token[1])));   //squareroot transform of intensity, try stabalized peak intentsity variance
				}
				line = bf.readLine();
			}
			if(line == null){
				return null;
			}
		}catch(IOException ioe){
			System.out.println("Cannot Open MGF file");
			System.out.println(ioe.getMessage());
			return null;
		}
		return s;
	}
	
	//need to fill in detail about parsing the annotation peptide information
	private static Spectrum readSpectrumFromMSP(BufferedReader bf) {
		Spectrum s = new SimpleSpectrum();
		try {
			String temp ;
			String line;
			int mod ;
	
			boolean isPeaks = false; 
			String[] tokens;
		   	line = bf.readLine() ;
		   
   			while (line != null && line.length() != 0) {
	   			
	   			if(line.startsWith("Name:")){
	   				s.setCharge(Integer.valueOf((line.split("[ ,/]"))[2]));
					//s.peptide = (line.split("[ ,/]"))[1];
					//s.peptide = s.peptide + "." + this.charge;
				
				}else if(line.startsWith("Comment:")){
					temp = line.substring(line.indexOf("Parent")) ;
					s.setParentMass(Double.parseDouble((temp.split("[=, ]"))[1])) ;
	
					temp = line.substring(line.indexOf("Mods")) ;
			    	tokens = temp.split("[,,=, ,/]") ;
			
					mod = Integer.parseInt(tokens[1]) ;
	    	
			    	if (mod == 0) {
			 
			    		//this.modPos = -1 ;
			    		//this.modMass = 0 ;
			    			    	
			    	}
			   
			    	else {
			    		
			    		//this.modPos = Integer.parseInt(tokens[2]) ;
			    		//this.modMass = 100.00 ;
			    			 
			    	}
			
				}else if(line.startsWith("Num peaks:")){
					isPeaks = true ;
				
				}
				
				else if(isPeaks) {
					tokens = line.split("\t");
					s.addPeak((new SimplePeak(Double.valueOf(tokens[0]),
							Double.valueOf(tokens[1]))));
				}
				line = bf.readLine();
			}
			
			if(line == null){
				return null;
			}
   		  			
   		}catch(IOException ioe){
			System.out.println("Cannot Open MST file");
			System.out.println(ioe.getMessage());
			return null;
		}
		return s;
	}
}
