package MSPLIT.SpectrumIO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.Spectrums.Peak;

import MSPLIT.SimpleAnnotatedSpectrum;
import MSPLIT.SimplePeak;
import MSPLIT.SimpleSpectrum;
import MSPLIT.Spectrum;


public class SpectrumReader {
	//MSP is one of the NIST spectral format so we should read in
	//a annotated spectrum
	public static Spectrum readSpectrumFromMSP(BufferedReader bf) {
		SimpleAnnotatedSpectrum s = new SimpleAnnotatedSpectrum();
		try {
		
			String temp ;
    		int size ;
			String line;
			int mod ;
	
			boolean isPeaks = false; 
			String[] tokens;
		   	line = bf.readLine() ;
		   
   			while (line != null && line.length() != 0) {
	   			if(line.startsWith("Name:")){
	   				s.setCharge(Integer.valueOf((line.split("[ ,/]"))[2]));
					s.setPeptide(line.split("[ ,/]")[1]);
					//this.peptide = this.peptide + "." + this.charge;
				
				}else if(line.startsWith("Comment:")){
					temp = line.substring(line.indexOf("Parent")) ;
	    			s.setParentMass(Double.parseDouble((temp.split("[=, ]"))[1])) ;
	
					temp = line.substring(line.indexOf("Mods")) ;
			    	tokens = temp.split("[,,=, ,/]") ;
			
					mod = Integer.parseInt(tokens[1]) ;
					//we will try to take care of PTM later
//			    	if (mod == 0) {			 
//			    		this.modPos = -1 ;
//			    		this.modMass = 0 ;
//			    			    	
//			    	}
//			   
//			    	else {
//			    		
//			    		this.modPos = Integer.parseInt(tokens[2]) ;
//			    		this.modMass = 100.00 ;
//			    			 
//			    	}
			
				}else if(line.startsWith("Num peaks:")){
					isPeaks = true ;
				
				}
				
				else if(isPeaks) {
					tokens = line.split("\t");
					s.addPeak(new SimplePeak(Double.valueOf(tokens[0]),
							Double.valueOf(tokens[1])));
				}
				line = bf.readLine();
				//System.out.println("line is: " + line);
			}
			
			if(line == null){
				//System.out.println("line is null");
				return null;
			}
			//System.out.println("returning spectrum");
   		  	return s;		
   		}catch(IOException ioe){
			System.out.println("Cannot Open MSP file");
			System.out.println(ioe.getMessage());
			return null;
		}
		
	}
	
	//again splib is a format of NIST spectral library should read in a 
	//annotated spectrum
	public static Spectrum readSpectrumFromSplib(BufferedReader bf) {
		SimpleAnnotatedSpectrum s = new SimpleAnnotatedSpectrum();
		try {
			String temp ;
    		int size ;
			String line;
			int mod ;
	
			boolean isPeaks = false; 
			String[] tokens;
		   	line = bf.readLine() ;
		   
   			while (line != null && line.length() != 0) {
	   			
	   			if(line.startsWith("Name:")){
	   				s.setCharge(Integer.valueOf((line.split("[ ,/]"))[2]));
					s.setPeptide(line.split("[ ,/]")[1]);
					//this.peptide = this.peptide + "." + this.charge;
					
				}else if(line.startsWith("PrecursorMZ:")){
	   				s.setParentMass(Double.parseDouble((line.split(" "))[1]));
	   				
				}else if(line.startsWith("Comment:")){	
//					if(line.contains("DECOY")){
//						this.peptide = "X_" + peptide;
//					}
//					temp = line.substring(line.indexOf("Mods")) ;
//			    	tokens = temp.split("[,,=, ,/]") ;
//			
//					mod = Integer.parseInt(tokens[1]) ;
//	    	
//			    	if (mod == 0) {
//			 
//			    		this.modPos = -1 ;
//			    		this.modMass = 0 ;
//			    			    	
//			    	}
//			   
//			    	else {
//			    		
//			    		this.modPos = Integer.parseInt(tokens[2]) ;
//			    		this.modMass = 100.00 ;
//			    			 
//			    	}
			
				}else if(line.startsWith("NumPeaks:")){
					isPeaks = true ;
				
				}
				
				else if(isPeaks) {
					tokens = line.split("\t");
					s.addPeak(new SimplePeak(Double.valueOf(tokens[0]),
							Double.valueOf(tokens[1])));
				}
				line = bf.readLine();
			}
			
			if(line == null){
				//System.out.println("line is null");
				return null;
			}
		System.out.println("read in peaks: " + s.numOfPeaks() + " for spectrum " + s.getPeptide());			
   		}catch(IOException ioe){
			System.out.println("Cannot Open splib file");
			System.out.println(ioe.getMessage());
			return null;
		}
   		//System.out.println("read in peaks: " + this.peaks.size() + " for spectrum " + this.peptide);
		return s;	
	}
	
	//ms2 is a simple format like the .dta, it do contain annotated information
	//for spectrum, probably should only create spectrums
	public static Spectrum readSpectrumFromMS2(BufferedReader bf) {
		SimpleSpectrum s = new SimpleSpectrum();
		try {
			String temp ;
    		int size ;
			String line;
			int mod ;
			boolean isPeaks = false; 
			String[] tokens;
		   	line = bf.readLine() ;
		   
   			while (line != null && !(isPeaks && line.startsWith("S"))) {
	   			
	   			if(line.startsWith("S")){
	   				tokens = line.split("\\s+");
	   				//this.scanNumber = Integer.parseInt(tokens[1]);
	   				//this.peptide = "Scan Number: " + this.scanNumber;
	   				s.setSpectrumName("Scan Number: " + Integer.parseInt(tokens[1]));
	   				s.setParentMass(Double.parseDouble(tokens[3]));
	   				s.setCharge(0); //unknown charge, contain no charge info???
				}else if(Character.isDigit(line.charAt(0))){
					isPeaks=true;
				}	
				if(isPeaks) {
					tokens = line.split("\\s+");
					s.addPeak(new SimplePeak(Double.valueOf(tokens[0]),
							Double.valueOf(tokens[1])));
				}
				bf.mark(200);
				line = bf.readLine();
			}
   			bf.reset(); //we backtrack the reader to previousline
			
			if(line == null){
				//System.out.println("line is null");
				return null;
			}
		System.out.println("read in peaks: " + s.numOfPeaks() + " for spectrum " + s.getSpectrumName());			
   		}catch(IOException ioe){
			System.out.println("Cannot Open splib file");
			System.out.println(ioe.getMessage());
			return null;
		}
   		//System.out.println("read in peaks: " + this.peaks.size() + " for spectrum " + this.peptide);
		return s;	
	}
	
	
	//read a spectrum from a buffer reader, this way
	//when we try to create multiple spectrums from MGF
	//we can pass the reader to it directly
	//note that means after calling this function we
	//read one subsequent spectra from the file
	//MGF format can contain both annotated/unannotated spectrum
	//so we will create annoated spectrum to account for all cases
	public static Spectrum readSpectrumFromMGF(BufferedReader bf){
		SimpleAnnotatedSpectrum s = new SimpleAnnotatedSpectrum();
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
					s.setParentMass(Double.valueOf(line.split("=")[1]));
					//isPeaks = true;  //for real data last line of mgf is pepmass
				}else if(line.startsWith("CHARGE")){
					charge = ((line.split("="))[1]);
					if(charge.startsWith("+")){
						charge = charge.substring(1);
					}
					s.setCharge(Integer.valueOf(charge));
					
				}else if(line.startsWith("PEPSEQ") || line.startsWith("SEQ")){
					s.setPeptide(line.split("=")[1]);
					//this.peptide = this.peptide + "." + this.charge;
				}else if(line.startsWith("TITLE")){
					String[] tokens = line.split("=");
					if(tokens.length > 1){
						//this.peptide = (line.split("="))[1];
						//this.peptide = this.peptide + "." + this.charge;
						s.setSpectrumName(line.split("=")[1]);
					}
				}else if(line.startsWith("PEPMOD")){
//					//just put some arbitary mod for now
//					//mainly use this as a flag to tell spectra with mod
//					//from those without any mod
//					this.modMass = 50;
//					this.modPos = 3;
					
				}else if(Character.isDigit(line.charAt(0))){
					//System.out.println("currelint is: " + line + "\tfirst: " + line.charAt(0));
					isPeaks = true;
					
				}
				if(isPeaks){
					token = line.split("\\s+");
					s.addPeak(new SimplePeak(Double.valueOf(token[0]), //note multiply by 0.9995 to make peaks center around integer values
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

}
