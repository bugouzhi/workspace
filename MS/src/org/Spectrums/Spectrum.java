//contain various informations about a spectrum 
package org.Spectrums;
import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Comparator;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Map.Entry;

public class Spectrum implements Comparable<Spectrum>, Serializable{
	//default value for vector representation
	public static final long serialVersionUID = 1L;
	private static double BINWIDTH = 1;
	private static double MINMASS = 0.5;
	private static double MAXMASS = 2000;
	//the peaks in this spectrum, we should always keep this list sorted
	//according to the mass
	private List <Peak> peaks;
	public String spectrumName = "";
	public String peptide;
	String protein;
	public double parentMass;
	public int charge;
	double modMass; //mass of any modification
	int modPos; //the position of modification
	double score = 0; //provide certain kind of score to this spectrum
	                  //we can utilize this field to order the spectrums
	double upperBound = 0;// upper bound for searching
	double rt = 0;
	public int scanNumber = 0; //scan number of the spectrum
	public Spectrum(){ //create dummy spectrum
			this.peaks =new Vector();
			this.spectrumName = "DUMMYSPEC";
			this.peptide = "DUMMYSEQ";
		    this.modMass = 0;
			this.modPos = 0;
			this.charge = 1;
			this.parentMass = 10000; //very large mass, so the dummy
			                         //is different than real spectra
			this.scanNumber = 0;
	}	
	
	public Spectrum(Spectrum s){
		this.peaks =new ArrayList<Peak>();
		this.spectrumName = s.spectrumName;
		this.peptide = s.peptide;
	    this.modMass = s.modMass;
		this.modPos = s.modPos;
		this.charge = s.charge;
		this.parentMass = s.parentMass;
		this.scanNumber = s.scanNumber;
		this.score = s.score;
		this.upperBound = s.upperBound;
		this.protein = s.protein;
		for(int i = 0; i < s.peaks.size(); i++){
			Peak p = s.peaks.get(i);
			//System.out.println("copying: " + p);
			this.peaks.add(new Peak(p));
			//System.out.println("copied: " + this.peaks.get(this.peaks.size()-1));
		}
		//Collections.sort(this.peaks, PeakMassComparator.comparator);
		//this.peaks = s.peaks;
	}
	//note in the peptide we also include the charge as well, since
	//we are using the peptide as key to group spectrum, thus spectrum
	//with same peptide AND same charge should be grouped together
	public Spectrum(Vector peaks, String peptide, double pm, int charge,
		double modMass, int modPos){
		this.peaks = peaks;
		this.peptide = peptide;
		this.parentMass = pm;
		this.charge = charge;
		this.modMass = modMass;
		this.modPos = modPos;
	}
	
	public Spectrum(Vector peaks, String peptide, String name, double pm, int charge,
			double modMass, int modPos, int scanNumber){
			this.peaks = peaks;
			this.peptide = peptide;
			this.parentMass = pm;
			this.charge = charge;
			this.modMass = modMass;
			this.modPos = modPos;
			this.spectrumName = name;
			this.scanNumber = scanNumber;
	}
	
	public Spectrum(String spectrumFile, String format){
		if(format.equals("MGF")){
			readSpectrumFromMGF(spectrumFile);
		}
	}
	
	public Spectrum readSpectrumFromNIST(String fileName){
		return new Spectrum();
	}
	
	//create a new spectrum by combining two spectrum together
	//note this is a simulated spectrum, so some of the field might
	//not make sense 
	public Spectrum(Spectrum s1, Spectrum s2, double scale1, double scale2){
		s1 = s1.toNormVector(BINWIDTH, MINMASS, MAXMASS);
		s2 = s2.toNormVector(BINWIDTH, MINMASS, MAXMASS);
		//s1 = s1.toVector(BINWIDTH, MINMASS, MAXMASS);
		//s2 = s2.toVector(BINWIDTH, MINMASS, MAXMASS);
		this.peptide = s1.peptide + " & " +  s2.peptide;
		this.modMass = 0;
		this.modPos = 0;
		this.charge = s1.charge;
		this.peaks = new Vector();
		this.parentMass = Math.max(s1.parentMass, s2.parentMass);  //we use the larger mass for the mixed spectrum
		//note here we make the modMass of a mixture the difference of pm of the components
		this.modMass = this.parentMass - Math.min(s1.parentMass, s2.parentMass);
		Peak p1, p2;
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < s1.peaks.size() && j < s2.peaks.size()){
			p1 = s1.peaks.get(i);
			p2 = s2.peaks.get(j);
			mz1 = p1.getMass();
			mz2 = p2.getMass();
			if(mz1 < mz2){
				this.peaks.add(new Peak(mz1, p1.getIntensity()*scale1));
				i++;
			}else if(mz1 == mz2){
				this.peaks.add(new Peak(mz1, p1.getIntensity()*scale1 + p2.getIntensity()*scale2));
				i++;
				j++;
			}else {
				this.peaks.add(new Peak(mz2, p2.getIntensity()*scale2));
				j++;
			}
		}
		//appending any remaining peaks 
		while(i < s1.peaks.size()){
			p1 = s1.peaks.get(i);
			this.peaks.add(new Peak(p1.getMass(), p1.getIntensity()*scale1));
			i++;
		}
		while(j < s2.peaks.size()){
			p2 = s2.peaks.get(j);
			this.peaks.add(new Peak(p2.getMass(), p2.getIntensity()*scale2));
			j++;
		}
	}
	
	//with copied rank and combine mass within tolerance
	public Spectrum(Spectrum s1, Spectrum s2, double scale1, double scale2, boolean toVector){
		if(toVector){
			s1 = s1.toVector(BINWIDTH, MINMASS, MAXMASS);
			s2 = s2.toVector(BINWIDTH, MINMASS, MAXMASS);
			s1.computePeakRank();
			s2.computePeakRank();
		}
		int commonpeak = 0;
		double mag1 = s1.sumOfPeaks();
		mag1 = mag1*mag1;
		double mag2 = s2.sumOfPeaks();
		mag2 = mag2*mag2;
		System.out.println(s1.peptide + " & " + s2.peptide + " real alpha " + (mag1/s1.magnitude())/(mag2/s2.magnitude()));
		//mag1 = s1.magnitude();
		//mag2 = s2.magnitude();
		s1.scaleSpectrum(1/mag1);
		s2.scaleSpectrum(1/mag2);
		if(s1.peptide.contains(".")){
			this.peptide = s1.peptide + " & " +  s2.peptide;
		}else{
			this.peptide = s1.peptide + "."+ s1.charge + " & " +  s2.peptide + "." + s2.charge;
		}
		this.modMass = 0;
		this.modPos = 0;
		this.charge = s1.charge + s2.charge;              //temporary change use to train model for mixture
		this.peaks = new Vector();
		this.parentMass = s1.parentMass;//Math.max(s1.parentMass, s2.parentMass);  //we use the larger mass for the mixed spectrum
		//note here we make the modMass of a mixture the difference of pm of the components
		this.modMass = this.parentMass - Math.min(s1.parentMass, s2.parentMass);
		Peak p1, p2;
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < s1.peaks.size() && j < s2.peaks.size()){
			p1 = s1.peaks.get(i);
			p2 = s2.peaks.get(j);
			mz1 = p1.getMass();
			mz2 = p2.getMass();
			if(mz1 < mz2){
				this.peaks.add(new Peak(mz1, p1.getIntensity()*scale1));
				peaks.get(peaks.size()-1).copyRank(p1);
				i++;
			}else if(mz1 == mz2){
				this.peaks.add(new Peak(mz1, p1.getIntensity()*scale1 + p2.getIntensity()*scale2));
				if(p1.getIntensity() > p2.getIntensity()){
					peaks.get(peaks.size()-1).copyRank(p1);
				}else{
					peaks.get(peaks.size()-1).copyRank(p2);
				}
				commonpeak++;
				i++;
				j++;
			}else {
				this.peaks.add(new Peak(mz2, p2.getIntensity()*scale2));
				peaks.get(peaks.size()-1).copyRank(p2);
				j++;
			}
		}
		//appending any remaining peaks 
		while(i < s1.peaks.size()){
			p1 = s1.peaks.get(i);
			this.peaks.add(new Peak(p1.getMass(), p1.getIntensity()*scale1));
			peaks.get(peaks.size()-1).copyRank(p1);
			i++;
		}
		while(j < s2.peaks.size()){
			p2 = s2.peaks.get(j);
			this.peaks.add(new Peak(p2.getMass(), p2.getIntensity()*scale2));
			peaks.get(peaks.size()-1).copyRank(p2);
			j++;
		}
		//rSystem.out.println("peaks count " + s1.getPeak().size() + "\t" + s2.getPeak().size() + "\t" + commonpeak);
	}
	
	public Spectrum(Spectrum s1, Spectrum s2, double scale1, double scale2, boolean toVector, boolean merge){
		this(s1, s2, scale1, scale2, toVector);
		if(merge){
			mergePeaks(this, 0.25);
		}
		
	}
	
	//iterativel merge peaks that are close enough in masses
	public void mergePeaks(Spectrum mix, double tolerance){
		//System.out.println("we start with  " + mix.peaks.size() + " peaks");
		
		for(int i = 0; i < mix.getPeaks().size(); i++){
			Peak current = mix.getPeaks().get(i);
			current.setIntensity(current.getIntensity()+ Math.random()*0.0000001);  //we add a very small number to make sure no intensity clashes
		}
		
		TreeMap<Double, Peak> intensityMap = new TreeMap<Double, Peak>();
		TreeMap<Double, Peak> massMap = new TreeMap<Double, Peak>();
		List<Peak> newPeakList = new ArrayList<Peak>();
		for(int i = 0; i < mix.getPeaks().size(); i++){
			Peak current = mix.getPeaks().get(i);
			current.setIntensity(current.getIntensity()+Math.random()*0.000000001);//in case exact same intensity
			current.setMoz(current.getMass()+Math.random()*0.000000001); //in case exact same mass
			intensityMap.put(current.getIntensity(), current); 
			massMap.put(current.getMass(), current); //in case exact same masses
		}
		while(intensityMap.size() > 0){
			Entry<Double, Peak> e = intensityMap.lastEntry();
			Peak current = e.getValue();
			Double currentKey = current.getIntensity();
			SortedMap<Double, Peak> subMap = massMap.subMap(current.getMass()-tolerance, current.getMass()+tolerance);
			//System.out.println("current: " + current);
			//System.out.println("map has size " + intensityMap.size());
			//System.out.println("map has size " + massMap.size());
			List<Peak> toBeRemoved = new ArrayList<Peak>();
			while(subMap.size() > 1){
				for(Iterator<Peak> it = subMap.values().iterator(); it.hasNext();){
					Peak neigh = it.next();
					if(neigh.equals(current)){
						continue;
					}
					double weight = current.getIntensity() / (current.getIntensity()+neigh.getIntensity());
					current.setMoz(current.getMass()*weight + neigh.getMass()*(1-weight));
					current.setIntensity(current.getIntensity()+neigh.getIntensity());
					toBeRemoved.add(neigh);
				}
				for(int j = 0; j < toBeRemoved.size(); j++){
					//System.out.println("add to removed: " + toBeRemoved.get(j));
					intensityMap.remove(toBeRemoved.get(j).getIntensity());
					massMap.remove(toBeRemoved.get(j).getMass());
					
				}
				subMap = massMap.subMap(current.getMass()-tolerance, current.getMass()+tolerance);
				//System.out.println("neighbor has size: " + subMap.size());
			}
			newPeakList.add(current);
			intensityMap.remove(currentKey);
			
			//System.out.println("remainging peaks to process: " + intensityMap.size());
		}
		Collections.sort(newPeakList, PeakMassComparator.comparator);
		mix.setPeaks(newPeakList);
		//System.out.println("after merging we have "  + mix.peaks.size());
	}
	
	public void deIsoPeaks(Spectrum mix, double tolerance){
		System.out.println("we start with  " + mix.peaks.size() + " peaks");
		for(int i = 0; i < mix.getPeaks().size(); i++){
			Peak current = mix.getPeaks().get(i);
			current.setIntensity(current.getIntensity()+ Math.random()*0.0000001);  //we add a very small number to make sure no intensity clashes
		}
		
		TreeMap<Double, Peak> intensityMap = new TreeMap<Double, Peak>();
		TreeMap<Double, Peak> massMap = new TreeMap<Double, Peak>();
		List<Peak> newPeakList = new ArrayList<Peak>();
		
		for(int i = 0; i < mix.getPeaks().size(); i++){
			Peak current = mix.getPeaks().get(i);
			current.setIntensity(current.getIntensity()+Math.random()*0.000000001);//in case exact same intensity
			current.setMoz(current.getMass()+Math.random()*0.000000001); //in case exact same mass
			intensityMap.put(current.getIntensity(), current); 
			massMap.put(current.getMass(), current); //in case exact same masses
		}
		
		double offset = Mass.C13 - Mass.C12;
		int numC13 = 2;
		while(intensityMap.size() > 0){
			Entry<Double, Peak> e = intensityMap.lastEntry();
			Peak current = e.getValue();
			Double currentKey = current.getIntensity();
			SortedMap<Double, Peak> subMap = massMap.subMap(current.getMass(), current.getMass()+tolerance + offset*numC13);
			//System.out.println("map has size " + intensityMap.size());
			//System.out.println("map has size " + massMap.size());
			List<Peak> toBeRemoved = new ArrayList<Peak>();
			if(subMap.size() > 1){
				for(Iterator<Peak> it = subMap.values().iterator(); it.hasNext();){
					Peak neigh = it.next();
					if(neigh.equals(current)){
						continue;
					}
					for(int c13 = 1; c13 <= numC13; c13++){
						if(Math.abs(neigh.getMass() - current.getMass() - offset) <tolerance){
							//current.setIntensity(current.getIntensity() + neigh.getIntensity());
							current.setRank(1);
							neigh.setRank(c13);
							toBeRemoved.add(neigh);
						}
					}
				}
				for(int j = 0; j < toBeRemoved.size(); j++){
					intensityMap.remove(toBeRemoved.get(j).getIntensity());
					massMap.remove(toBeRemoved.get(j).getMass());
				}
				//subMap = massMap.subMap(current.getMass()-tolerance, current.getMass()+tolerance);
				//System.out.println("neighbor has size: " + subMap.size());
			}
			newPeakList.add(current);
			intensityMap.remove(currentKey);
			//System.out.println("remainging peaks to process: " + intensityMap.size());
		}
		Collections.sort(newPeakList, PeakMassComparator.comparator);
		mix.setPeaks(newPeakList);
		//System.out.println("after merging we have "  + mix.peaks.size());
	}
	
	
	
	
	
	public double getMeanIntensity(List<Peak> pList){
		double sum = 0.0;
		for(int i = 0; i < pList.size(); i++){
			sum += pList.get(i).getIntensity();
		}
		return sum/pList.size();
	}
	
	public double getIntensityDeviation(List<Peak> pList, double mean){
		double variance = 0.0;
		double deviation = 0.0;
		for(int i = 0; i < pList.size(); i++){
			deviation = pList.get(i).getIntensity() - mean;
			variance += deviation*deviation;
		}
		return Math.pow(variance, 0.5)/pList.size();
	}
	
	
	//same as above except this time we are not scaling the intensity
	public Spectrum(Spectrum s1, Spectrum s2){
		this(s1,s2,1,1);
	}
	
	public void setPeaks(List<Peak> peaks){
		this.peaks = peaks;
	}
	
	public List<Peak> getPeak(){
		return this.peaks;
	}
	
	public void readSpectrumFromMSP(String fileName){
		try{
			BufferedReader bf = new BufferedReader(new FileReader(fileName));
			readSpectrumFromMSP(bf);
		 }catch(IOException ioe){
				System.out.println("Cannot Open MSP file");
				System.out.println(ioe.getMessage());			
		 }
	}
	
	public boolean readSpectrumFromMSP(BufferedReader bf) {

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
	   				this.charge = Integer.valueOf((line.split("[ ,/]"))[2]);
					this.peptide = (line.split("[ ,/]"))[1];
				
				}else if(line.startsWith("Comment:")){
					temp = line.substring(line.indexOf("Parent")) ;
	    			this.parentMass = Double.parseDouble((temp.split("[=, ]"))[1] ) ;
	
					temp = line.substring(line.indexOf("Mods")) ;
			    	tokens = temp.split("[,,=, ,/]") ;
			
					mod = Integer.parseInt(tokens[1]) ;
	    	
			    	if (mod == 0) {
			 
			    		this.modPos = -1 ;
			    		this.modMass = 0 ;
			    			    	
			    	}
			   
			    	else {
			    		
			    		this.modPos = Integer.parseInt(tokens[2]) ;
			    		this.modMass = 100.00 ;
			    			 
			    	}
			
				}else if(line.startsWith("Num peaks:")){
					isPeaks = true ;
				
				}
				
				else if(isPeaks) {
					tokens = line.split("\t");
					this.peaks.add(new Peak(Double.valueOf(tokens[0]),
							Double.valueOf(tokens[1])));
				}
				line = bf.readLine();
			}
			
			if(line == null){
				//System.out.println("line is null");
				return false;
			}
   		  			
   		}catch(IOException ioe){
			System.out.println("Cannot Open MSP file");
			System.out.println(ioe.getMessage());
			return false;
		}
		return true;
		
		
	}
	
	public boolean readSpectrumFromSplib(BufferedReader bf) {
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
	   				this.charge = Integer.valueOf((line.split("[ ,/]"))[2]);
					String pep = (line.split("[ ,/]"))[1];
	   				this.peptide = new String(pep);
					//this.peptide = this.peptide;// + "." + this.charge;
				
				}else if(line.startsWith("PrecursorMZ:")){
	   				this.parentMass = Double.parseDouble((line.split(" "))[1]);
					
				}else if(line.startsWith("Comment:")){
					temp = line.substring(line.indexOf("Protein")) ;
			    	String prot = temp.split("[=\\s+]")[1] ;
			    	if(prot.startsWith("\"")){  //remove quotation at the beginning of protein name
			    		prot = prot.substring(1);
			    	}
					this.protein = new String(prot);		    	
					if(line.contains("DECOY")){
						//this.peptide = "X_" + peptide;
						this.protein = "DECOY_"+this.protein;
					}
					temp = line.substring(line.indexOf("Mods")) ;
			    	tokens = temp.split("[,,=, ,/]") ;
			
					mod = Integer.parseInt(tokens[1]) ;
	    	
			    	if (mod == 0) {			 
			    		this.modPos = -1 ;
			    		this.modMass = 0 ;
			    	}
			   
			    	else {
			    		
			    		this.modPos = Integer.parseInt(tokens[2]) ;
			    		this.modMass = 100.00 ;
			    			 
			    	}
			
				}else if(line.startsWith("NumPeaks:")){
					int numPeaks = Integer.parseInt(line.split(": ")[1]);
					isPeaks = true ;
					this.peaks = new ArrayList(numPeaks);
				}else if(isPeaks) {
					tokens = line.split("\t");
					this.peaks.add(new Peak(Double.valueOf(tokens[0]),
							Double.valueOf(tokens[1])));
				}
				line = bf.readLine();
			}
			if(line == null){
				//System.out.println("line is null");
				return false;
			}
		//System.out.println("read in peaks: " + this.peaks.size() + " for spectrum " + this.peptide);			
   		}catch(IOException ioe){
			System.out.println("Cannot Open splib file");
			System.out.println(ioe.getMessage());
			return false;
		}
   		//System.out.println("read in peaks: " + this.peaks.size() + " for spectrum " + this.peptide);
		return true;	
	}
	
	public boolean readSpectrumFromMS2(BufferedReader bf) {
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
	   				this.scanNumber = Integer.parseInt(tokens[1]);
	   				this.peptide = "Scan Number: " + this.scanNumber;
	   				this.spectrumName = "Scan Number: " + this.scanNumber;
	   				this.parentMass = Double.parseDouble(tokens[3]);
	   				this.charge = 0; //unknown charge
				}else if(Character.isDigit(line.charAt(0))){
					isPeaks=true;
				}	
				if(isPeaks) {
					tokens = line.split("\\s+");
					this.peaks.add(new Peak(Double.valueOf(tokens[0]),
							Double.valueOf(tokens[1])));
				}
				bf.mark(200);
				line = bf.readLine();
			}
   			bf.reset(); //we backtrack the reader to previousline
			
			if(line == null){
				//System.out.println("line is null");
				return false;
			}
		System.out.println("read in peaks: " + this.peaks.size() + " for spectrum " + this.peptide);			
   		}catch(IOException ioe){
			System.out.println("Cannot Open splib file");
			System.out.println(ioe.getMessage());
			return false;
		}
   		//System.out.println("read in peaks: " + this.peaks.size() + " for spectrum " + this.peptide);
		return true;	
	}
	
	//read a spectrum from a MGF file
	public void readSpectrumFromMGF(String fileName){
		try{
			BufferedReader bf = new BufferedReader(new FileReader(fileName));
			readSpectrumFromMGF(bf);
		 }catch(IOException ioe){
				System.out.println("Cannot Open MGF file");
				System.out.println(ioe.getMessage());			
		 }
	}
	
	//read a spectrum from a buffer reader, this way
	//when we try to create multiple spectrums from MGF
	//we can pass the reader to it directly
	//note that means after calling this function we
	//read one subsequent spectra from the file
	public boolean readSpectrumFromMGF(BufferedReader bf){
		try{
			String line;
			//System.out.println("begin parsing");
			do{
				line = bf.readLine();
				//System.out.println("line is: " + line);
			}while(line != null && !line.contains("BEGIN IONS"));
			// a flag tell us whether we are in peak sections
            //in the spectrum file
			boolean isPeaks = false; 
			String[] token;
			String charge;
			line = bf.readLine();
			while(line != null && !line.contains("END IONS")){
				//System.out.println("line is " + line);
				if(line.startsWith("PEPMASS")){
					String[] mass = line.split("=");
					mass = mass[1].split("\\s+");
					this.parentMass = Double.valueOf(mass[0]);
					//isPeaks = true;  //for real data last line of mgf is pepmass
				}else if(line.startsWith("CHARGE")){
					charge = ((line.split("[=\\s*]"))[1]);
					if(charge.startsWith("+")){
						charge = charge.substring(1);
					}
					if(charge.contains("+")){
						charge=charge.replace("+", "");
					}
					this.charge = Integer.valueOf(charge);
					
				}else if(line.startsWith("PEPSEQ") || line.startsWith("SEQ")){
					this.peptide = (line.split("[=]"))[1];
					//this.peptide = this.peptide + "." + this.charge;
				}else if(line.startsWith("SCAN")){
					this.scanNumber = (Integer.parseInt(line.split("[=]")[1]));
					//this.peptide = this.peptide + "." + this.charge;
				}else if(line.startsWith("TITLE")){
					String[] tokens = line.split("=");
//					if(tokens.length > 1){
//						this.peptide = (line.split("="))[1];
//						//this.peptide = this.peptide + "." + this.charge;
//						this.spectrumName = (line.split("="))[1];
//						if(this.spectrumName.contains("Scan Number:")){
//							String[] tokens2 = this.spectrumName.split("\\s+");
//							this.scanNumber = Integer.parseInt(tokens2[2]);
//						}
//					}
					this.spectrumName = line.substring(6);
				}else if(line.startsWith("PEPMOD")){
					//just put some arbitary mod for now
					//mainly use this as a flag to tell spectra with mod
					//from those without any mod
					this.modMass = 0;
					this.modPos = 0;
					
				}else if(Character.isDigit(line.charAt(0))){
					//System.out.println("currelint is: " + line + "\tfirst: " + line.charAt(0));
					isPeaks = true;
					
				}
				if(isPeaks){
					token = line.split("\\s+");
					//System.out.println("line is: " + line);
					try{
						//if(token.length < 2){
						//	System.out.println("line is: " + line);
						//}
						this.peaks.add(new Peak(Double.valueOf(token[0]), //note multiply by 0.9995 to make peaks center around integer values
						Double.valueOf(token[1])));   //
					}catch(NumberFormatException e){
						break;//return whatever is readed
					}
					
				}
				line = bf.readLine();
			}
			List<Peak> pListSlim = new ArrayList(this.peaks.size());
			pListSlim.addAll(this.peaks);
			this.peaks = pListSlim;
			if(line == null){
				return false;
			}
		}catch(IOException ioe){
			System.out.println("Cannot Open MGF file");
			System.out.println(ioe.getMessage());
			return false;
		}
		return true;
	}
	
	public int compareTo(Spectrum s){
		if(s.score == this.score){
			return 0;
		}else if(s.score < this.score){
			return 1;
		}else{
			return -1;
		}
	}
	//return a vector representation of the spectrum
	//where we bin the spectrums and sum up all the 
	//operationally this is just another spectrum with
	//whose x-coordinate of peaks are in the unit of bin-width
	//rather than in dalton, so it each peak will have 
	//an integer as their x-coordinate, e.g (1, 100.2) that means
	//we have a sum of peaks with total intensity 100.2 at bin 1 etc..
	public Spectrum toNormVector(double binWidth, double minMass, double maxMass){
		Spectrum newSpectrum = toVector(binWidth, minMass, maxMass);
		newSpectrum.normalize();
		return newSpectrum;
	}
	
	public Spectrum toNormVector(){
		
		return this.toNormVector(BINWIDTH, MINMASS, MAXMASS);
	}
	
	public Spectrum toVector(double binWidth, double minMass, double maxMass){
		int bins = (int)((maxMass-minMass)/binWidth) + 1;
		Vector <Peak>  newPeaks = new Vector();
		double rightBoundary = minMass + binWidth;
		double currentValue;
		double moz;
		int j = 0;
		int i = 0;
		//System.out.println("bins: " + bins);
  		for(i = 0; i < bins; i++){
			rightBoundary = minMass + i*binWidth;
			currentValue = 0;
			while(j < this.peaks.size()){
				moz = ((Peak)peaks.get(j)).getMass();
				//System.out.println("moz " + moz);
				//System.out.println("r-edge: " + rightBoundary);
				if(moz > rightBoundary){
					break;
				}else{
					currentValue += peaks.get(j).getIntensity(); 
					j++;
				}	
			}
			if(currentValue > 0){
				//System.out.println("creating new peaks");
				newPeaks.add(new Peak(i, currentValue));
			}
		}
		return new Spectrum(newPeaks, this.peptide, this.parentMass, this.charge,
				this.modMass, this.modPos);
		
	}
	
	//we normalize the intensity of the spectrum
	//so it is has norm one
	private void normalize(){
		double total = magnitude();
		this.scaleSpectrum(1/total);
	}
	
	//convert the spectrum to relative intensity
	public void toRelIntensity(){
		double biggest = getMaxIntensity(); 
		this.scaleSpectrum(1/biggest);
	}
	
	public String getPeptide() {
		return peptide;
	}
	public void setPeptide(String peptide) {
		this.peptide = peptide;
	}
	public List<Peak> getPeaks() {
		return peaks;
	}
	public double getMaxIntensity(){
		double max = 0.0;
		for(int i = 0; i < this.peaks.size(); i++){
			if(peaks.get(i).getIntensity() > max){
				max = peaks.get(i).getIntensity();
			}
		}
		return max;
	}
	
	//if a intensity gets two large 
	//e.g. >90% of total intensity we 
	public boolean checkOutLiner(){
		int i = 0;
		double total = magnitude();
		for(i = 0; i < this.peaks.size(); i++){
			if(peaks.get(i).getIntensity() > 0.9*total){
				return true;
			}
		}
		return false;
	}
	//return the magnitue of this spectrum
	//now we are thinking a spectrum as a vector
	public double magnitude(){
		double total = 0;
		for(int i = 0; i < peaks.size(); i++){
			total += (peaks.get(i)).getIntensity() * (peaks.get(i)).getIntensity();
		}
		total = Math.pow(total, 0.5);
		return total;
	}
	
	
	//magnitude for sqrt the spectrum 
	public double sumOfPeaks(){
		double total = 0;
		for(int i = 0; i < peaks.size(); i++){
			total += (peaks.get(i)).getIntensity();
		}
		total = Math.pow(total, 0.5);
		return total;
	}
	
	//magnitude by summing all peak intensity 
	public double sumMagnitude(){
		double total = 0;
		for(int i = 0; i < peaks.size(); i++){
			total += (peaks.get(i)).getIntensity();
		}
		return total;
	}
	
	//we scale all the peaks in this spectrum by some factor
	public void scaleSpectrum(double factor){
		for(int i = 0; i < peaks.size(); i++){
			peaks.get(i).scaleIntensity(factor);
		}
	}
	
	public void scaleMass(double factor){
		for(int i = 0; i < peaks.size(); i++){
			peaks.get(i).scaleMass(factor);
		}
	}
	
	//shift all masses in the spectrum by ceratin small mass
	public void shiftSpectrum(double shift){
		for(int i = 0; i < peaks.size(); i++){
			peaks.get(i).shiftMass(shift);
		}
	}
	
	public void shiftSpectrumPPM(double shiftPPM){
		for(int i = 0; i < peaks.size(); i++){
			peaks.get(i).shiftMassPPM(shiftPPM);
		}
		this.parentMass = this.parentMass + shiftPPM*parentMass/1000000;
	}
	
	//taking square root of all the intensity in spectrum
	//try to dampen the dominant effect of very high intensity
	public void sqrtSpectrum(){
		for(int i = 0; i < peaks.size(); i++){
			peaks.get(i).setIntensity(Math.pow(peaks.get(i).getIntensity(), 0.5));
		}
	}
	//compare two spectrum by calculating their normalized dot product 
	//i.e. cosine of vector representation of spectrum
	public double cosineSim(Spectrum s1){
		double product = 0;
		double magnitude = this.magnitude(); 
		magnitude *= s1.magnitude();
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += this.peaks.get(i).getIntensity()
					* s1.peaks.get(j).getIntensity();
				i++;
				j++;
			}else{
				j++;
			}
		}
		
		return product/magnitude;
		
	}
	
	public double cosine(Spectrum s1){
		return 0.0;
	}
	
	//for determine whether spectral similarity is dominate by
	//a few dominant peaks or many noise peaks
	public double cosineSim2(Spectrum s1){
		double product = 0;
		double magnitude = this.magnitude(); 
		//System.out.println("magnitue one: " + magnitude);
		//System.out.println("magnitue two: " + s1.magnitude());
		magnitude *= s1.magnitude();
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			//System.out.println(this.peaks.get(i));
			//System.out.println(s1.peaks.get(j));
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += Math.pow(this.peaks.get(i).getIntensity(), 2)
					* Math.pow(s1.peaks.get(j).getIntensity(), 2);
				i++;
				j++;
			}else{
				j++;
			}
		}
		return Math.pow(product, 0.5)/magnitude;	
	}
	
	//consine sim with sqrt normalization of peak intensities
	public double cosineSim1(Spectrum s1){
		double product = 0;
		double magnitude = this.sumOfPeaks(); 
		//System.out.println("magnitue one: " + magnitude);
		//System.out.println("magnitue two: " + s1.sumOfPeaks());
		magnitude *= s1.sumOfPeaks();
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += Math.pow(this.peaks.get(i).getIntensity(), 0.5)
					* Math.pow(s1.peaks.get(j).getIntensity(), 0.5);
				i++;
				j++;
			}else{
				j++;
			}
		}
		
		return product/magnitude;
		
	}
	
	//this method consider similarity of two spectrum
	//and their shifted variants, as different isotopic 
	//variants of same spectrum might have mass off 
	//by a small amount of mass/charge 
	//note this version should apply to the raw spectrum
	//i.e. without any normalization
	public double shiftCosineSim(Spectrum s1){
		Spectrum original = this.toNormVector();
		this.shiftSpectrum(0.5);
		Spectrum rightShift = this.toNormVector();
		this.shiftSpectrum(-1.0); //left shift original spectrum by -0.5 (0 + 0.5 - 1.0) = -0.5
		Spectrum leftShift = this.toNormVector(); 
		this.shiftSpectrum(0.5); //restore the original spectrum
		Spectrum vS1 = s1.toNormVector();
		double cosineOriginal = original.cosineSim(vS1);
		double cosineLeftShift = leftShift.cosineSim(vS1);
		double cosineRightShift = rightShift.cosineSim(vS1);
		//System.out.println("three version of cosine: " + cosineOriginal + ", " + cosineLeftShift + ", " + cosineRightShift);
		return Math.max(cosineRightShift, 
			Math.max(cosineOriginal, cosineLeftShift));
		
	}
	
	public double shareSim(Spectrum s1){
		double s = this.shareSim(s1, 0);
		return s;
	}
	
	//we calculate fraction of peaks share by the two spectrum,
	//there are three way we can calculate the fraction, which maybe useful at diff situation
	//mode: 0 , divide by average size of two spectrum being compare
	//mode: 1, divide by size of this spectrum
	//mode: 1, divide by size of s1
	public double shareSim(Spectrum s1, int mode){
		double product = 0;
		double magnitude = this.peaks.size(); 
		double magnitude2 = s1.peaks.size();
		//System.out.println("mixture magnitue " + s1.magnitude());
		double mz1, mz2; 
		int i = 0, j = 0;
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				product += 1;
				i++;
				j++;
			}else{
				j++;
			}
		}
		if(mode == 0){
			return product/((magnitude+magnitude2)/2);
		}else if(mode == 1){
			return product/magnitude;
		}else{
			return product/magnitude2;
		}
		
	}
	
	//similar to calculating cosine of two vector
	//but here we consider only those bins
	//that are non-zero in this vector
	public double projectedCosine(Spectrum s1){
		double product = 0;
		double magnitude = this.magnitude();
		//we do the dotproduct as above, but just calculates
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector
		double projectedNorm = 0.00000001; //very small number avoid div-by-zero error  
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += this.peaks.get(i).getIntensity()
					* s1.peaks.get(j).getIntensity();
				projectedNorm += s1.peaks.get(j).getIntensity() * s1.peaks.get(j).getIntensity();
				i++;
				j++;
			}else{
				j++;
			}
		}
		//System.out.println("this norm: " + this.magnitude());
		//System.out.println("mixture projected norm: " + Math.pow(projectedNorm, 0.5));
		magnitude = magnitude * Math.pow(projectedNorm, 0.5);
		return product/magnitude;
	}
	
	
	public double cosine(Spectrum s1, double tolerance){
		double product = 0;
		double magnitude = this.magnitude();
		double magnitude2 = s1.magnitude();
		//we do the dotproduct as above, but just calculates
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			//System.out.println(i+ "\t" +j);
			if(mz2 - mz1 > tolerance){
				//System.out.println("missing " + mz1 + "\t" + this.peaks.get(i).getIntensity() +'\t' + 0 + "\t" + 0);
				i++;
				continue;
			}
				
			if(mz1 - mz2 > tolerance){
				//System.out.println("missing " +  0 + "\t" + 0 + "\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity());
				j++;
				continue;
			}
			
			int p = i;
			Peak best=null;
			double bestIntensity = -1000;
			while(p < this.peaks.size()){
				double currMass = this.peaks.get(p).getMass();
				double currInt = this.peaks.get(p).getIntensity();
				if(currMass - mz2 > tolerance){
				//if(!SWATHUtils.checkMass(currMass, mz2, 32, SWATHUtils.PPM)){
					break;
				}
				//System.out.println("matched " + currMass + "\t" + this.peaks.get(p).getIntensity() + "\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity());
				//bestIntensity = bestIntensity < currInt ? currInt : bestIntensity;
				product += this.peaks.get(p).getIntensity()
					* s1.peaks.get(j).getIntensity();
				p++;
			}
			j++;
		}
	
		magnitude = magnitude * magnitude2;		
		//System.out.println("product: " + product  + "\t"  + "proj-norm: " + magnitude);
		return product/magnitude;
	}
	
	//an approximation of cosine using arbitrary tolerance
	//it is not exactly cosine since the magnitude computation is technically not exact
	public double cosineApprox(Spectrum s1, double tolerance){
		double product = 0;
		double magnitude = this.magnitude();
		double magnitude2 = s1.magnitude();
		//we do the dotproduct as above, but just calculates
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			//System.out.println(i+ "\t" +j);
			if(mz2 - mz1 > tolerance){
				i++;
				continue;
			}
			
			if(mz1 - mz2 > tolerance){
				j++;
				continue;
			}
			
			int p = i;
			Peak best=null;
			double bestIntensity = -1000;
			while(p < this.peaks.size()){
				double currMass = this.peaks.get(p).getMass();
				double currInt = this.peaks.get(p).getIntensity();
				if(currMass - mz2 > tolerance){
					break;
				}
				//System.out.println("matched " + currMass + "\t" + this.peaks.get(p).getIntensity() + "\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity());
				//bestIntensity = bestIntensity < currInt ? currInt : bestIntensity;
				product += this.peaks.get(p).getIntensity()
				* s1.peaks.get(j).getIntensity();
				p++;
			}
			//product += bestIntensity
			//* s1.peaks.get(j).getIntensity();
			
			j++;
		}
		//System.out.println("matched:" + product + "\tmagnitude:\t" + magnitude*magnitude2);
		//System.out.println("mixture projected norm: " + Math.pow(projectedNorm, 0.5));
		magnitude = magnitude * magnitude2;
		return product/magnitude;
	}
	
	public double projectedCosine(Spectrum s1, double tolerance){
		double product = 0;
		double magnitude = this.magnitude();
		//we do the dotproduct as above, but just calculates
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector
		double projectedNorm = 0.00000001; //very small number avoid div-by-zero error  
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			//System.out.println(i+ "\t" +j);
			if(mz2 - mz1 > tolerance){
				//System.out.println("missing " + mz1 + "\t" + this.peaks.get(i).getIntensity() +'\t' + 0 + "\t" + 0);
				i++;
				continue;
			}
				
			if(mz1 - mz2 > tolerance){
				j++;
				continue;
			}
			
			int p = i;
			Peak best=null;
			double bestIntensity = -1000;
			while(p < this.peaks.size()){
				double currMass = this.peaks.get(p).getMass();
				double currInt = this.peaks.get(p).getIntensity();
				//if(currMass - mz2 > tolerance){
				if(!SWATHUtils.checkMass(currMass, mz2, 46, SWATHUtils.PPM)){
					break;
				}
				//System.out.println("matched " + currMass + "\t" + this.peaks.get(p).getIntensity() + "\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity());
				//bestIntensity = bestIntensity < currInt ? currInt : bestIntensity;
				product += this.peaks.get(p).getIntensity()
					* s1.peaks.get(j).getIntensity();
				projectedNorm += s1.peaks.get(j).getIntensity() * s1.peaks.get(j).getIntensity();
				p++;
			}
			j++;
		}
		//System.out.println("this norm: " + this.magnitude());
	
		magnitude = magnitude * Math.pow(projectedNorm, 0.5);
		//System.out.println("product: " + product  + "\t"  + "proj-norm: " + magnitude);
		return product/magnitude;
	}
	
	//we are allowed to skip some peaks if there are intereference
	public double projectedCosineWithSkip(Spectrum s1, double tolerance, int maxGap){
		double pcos1 = this.projectedCosine(s1, tolerance);
		Spectrum projSpect = s1.project(this, tolerance);
		Spectrum copy = new Spectrum(this);
		copy.normalize();
		projSpect.normalize();
		double maxDiff = 0;
		Peak interference = null;
		Iterator<Peak> it1 = copy.peaks.iterator();
		Iterator<Peak> it2 = projSpect.peaks.iterator();
		if(!it2.hasNext()){
			return pcos1;
		}
		Peak p = it1.next();
		Peak p2 = it2.next();
		while(it1.hasNext() && it2.hasNext()){
			if(Math.abs(p.getMass() - p2.getMass()) < tolerance){
				double diff = Math.abs(p.getIntensity() - p2.getIntensity());
				//System.out.println("diff: " + diff + "\t" + p.getMass());
				if(diff > maxDiff){
					maxDiff = diff;
					interference = p;
				}
			}
			if(p.getMass() < p2.getMass()){
				p=it1.next();
			}else{
				p2=it2.next();
			}
		}
		//if(interference != null)	
		//	System.out.println("skipping: " + interference.getMass());
		copy.peaks.remove(interference);
		double pcos2 = copy.projectedCosine(s1, tolerance);
		//System.out.println("comparing: " + this.peptide + s1.peptide);
		//System.out.println("Detecting intereference using skip1 proj-cosine: " + pcos2 + "\tdelta:\t" + (pcos2-pcos1));
		if(pcos2 > pcos1){
			return pcos2;
		}else{
			return pcos1;
		}
	}
	
	
	//old-version....match first peak within tolerance... not neccesarilly correct 
	//but should be almost correct with merging of peaks within tolerance in spectrum
	public double projectedCosine2(Spectrum s1, double tolerance){
		double product = 0;
		double magnitude = this.magnitude();
		//we do the dotproduct as above, but just calculates
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector
		double projectedNorm = 0.00000001; //very small number avoid div-by-zero error  
		double mz1, mz2; 
		int i = 0, j = 0;
		System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 - mz2 < -1*tolerance){
				i++;
			}else if(mz1 - mz2 < tolerance){
				System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity());
				product += this.peaks.get(i).getIntensity()
					* s1.peaks.get(j).getIntensity();
				projectedNorm += s1.peaks.get(j).getIntensity() * s1.peaks.get(j).getIntensity();
				i++;
				j++;
			}else{
				j++;
			}
		}
		//System.out.println("this norm: " + this.magnitude());
		//System.out.println("mixture projected norm: " + Math.pow(projectedNorm, 0.5));
		magnitude = magnitude * Math.pow(projectedNorm, 0.5);
		return product/magnitude;
	}
	
	
	public double sharePeaks(Spectrum s1, double tolerance){
		return sharePeaks(s1, tolerance, false);
	}
	
	//share peak counts
	public double sharePeaks(Spectrum s1, double tolerance, boolean detail){
		double shareCount = 0;
		//we do the dotproduct as above, but just calculates
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector
		double mz1, mz2; 
		int i = 0, j = 0;
		if(detail)
			System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 - mz2 < -1*tolerance){
				if(detail)
					System.out.println("missing " + mz1 + "\t" + this.peaks.get(i).getIntensity() +"\t" + 0.00000);
				i++;
			}else if(mz1 - mz2 < tolerance){
				if(detail)
					System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				shareCount++;
				i++;
				j++;
			}else{
				j++;
			}
		}
		//System.out.println("this norm: " + this.magnitude());
		//System.out.println("mixture projected norm: " + Math.pow(projectedNorm, 0.5));
		//return shareCount/magnitude;
		return shareCount;
	}
	

	
	//compute number of shared peaks and number of shared peaks that are consider
	//signal, signal peaks are defined as a statistic in relation to noise intensity
	//currently it is implemented as ratio of intensity to mean intensity of noises
	public double[] projectedShare(Spectrum s1, double tolerance, double minZ){
		int shareCount=0;
		int bigCount = 0;
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		boolean isMatched = false;
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			//System.out.println(i+ "\t" +j);
			if(mz2 - mz1 > tolerance){
				i++;
				continue;
			}
			
			if(mz1 - mz2 > tolerance){
				if(!isMatched){
					System.out.println(s1.peptide + "\tmissing\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity() + "\t" + s1.peaks.get(j).getRank());
				}else{
					isMatched = false;
				}
				j++;
				
				continue;
			}
			
			int p = i;
			Peak best=null;
			double bestIntensity = -1000;
			boolean localShare = false;
			boolean localBig = false;
			while(p < this.peaks.size()){
				double currMass = this.peaks.get(p).getMass();
				double currInt = this.peaks.get(p).getIntensity();
				if(currMass - mz2 > tolerance){
					break;
				}
				System.out.println(s1.peptide + "\tmatched\t" + currMass + "\t" + this.peaks.get(p).getIntensity() + "\t" + this.peaks.get(p).getRank() +"\t"
						+ mz2 + "\t" + s1.peaks.get(j).getIntensity() + "\t" + s1.peaks.get(j).getRank());
				//bestIntensity = bestIntensity < currInt ? currInt : bestIntensity;
				localShare=true;
				isMatched=true;
				if(this.peaks.get(p).getRank() >= minZ){
				//if(this.peaks.get(p).getIntensity() > minZ){
					localBig=true;
				}
				p++;
			}
			if(localShare){
				shareCount++;
			}
			
			if(localBig) bigCount++;
			j++;
		}
		
		return new double[]{shareCount, bigCount};
	}
	
	//compute consensus of two spectrum, if peaks appear in both
	//increase the peak ranks, here rank is use as place holder
	//to represent the number of neighboring scans they appear
	public void computeConsensus(Spectrum s1, double tolerance){
		int shareCount=0;
		int bigCount = 0;
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			//System.out.println(i+ "\t" +j);
			if(mz2 - mz1 > tolerance){
				i++;
				continue;
			}
			
			if(mz1 - mz2 > tolerance){
				j++;
				continue;
			}
			
			int p = i;
			Peak best=null;
			double combinedIntensity = 0;
			boolean localShare = false;
			boolean localBig = false;
			while(p < this.peaks.size()){
				double currMass = this.peaks.get(p).getMass();
				double currInt = this.peaks.get(p).getIntensity();
				if(currMass - mz2 > tolerance){
					break;
				}
				//System.out.println("matched " + currMass + "\t" + this.peaks.get(p).getIntensity() + "\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity());
				//bestIntensity = bestIntensity < currInt ? currInt : bestIntensity;
				combinedIntensity += this.peaks.get(p).getIntensity();
				localShare=true;
				p++;
			}
			if(localShare){
				s1.peaks.get(j).setRank(s1.peaks.get(j).getRank()+1);
				//s1.peaks.get(j).setIntensity(s1.peaks.get(j).getIntensity() + combinedIntensity);
			}
			j++;
		}
	}
	
	//the sum of intensity of the mathed peaks
	//if there are more than one match can use to compute relative abundance
	public double projectedPeakIntensity(Spectrum s1, double tolerance){
		double shareIntensity = 0;
		double magnitude = s1.magnitude();
		magnitude=magnitude*magnitude;
		//we do the dotproduct as above, but just calculates
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector		
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			//System.out.println(i+ "\t" +j);
			if(mz2 - mz1 > tolerance){
				i++;
				continue;
			}
			
			if(mz1 - mz2 > tolerance){
				j++;
				continue;
			}
			
			int p = i;
			Peak best=null;
			double bestIntensity = -1000;
			while(p < this.peaks.size()){
				double currMass = this.peaks.get(p).getMass();
				double currInt = this.peaks.get(p).getIntensity();
				if(currMass - mz2 > tolerance){
					break;
				}
				//System.out.println("matched " + currMass + "\t" + this.peaks.get(p).getIntensity() + "\t" + mz2 + "\t" + s1.peaks.get(j).getIntensity());
				//bestIntensity = bestIntensity < currInt ? currInt : bestIntensity;
				shareIntensity += s1.peaks.get(j).getIntensity()*s1.peaks.get(j).getIntensity();
				break;
				//p++;
			}
			j++;
		}
//		System.out.println("this norm: " + magnitude);
		//System.out.println("mixture projected norm: " + Math.pow(projectedNorm, 0.5));
		//magnitude = magnitude * Math.pow(projectedNorm, 0.5);
		return shareIntensity/magnitude;
	}
	
	
	//similar to calculating cosine of two vector
	//but here we consider only those bins
	//that are non-zero in this vector
	public double projectedCosine1(Spectrum s1){
		double product = 0;
		double magnitude = this.sumOfPeaks();
		//we do the dotproduct as above, but just calculate
		//the norm for s1 vector differently, we only consider
	 	//those values that are non-zero in this vector
		double projectedNorm = 0.00000001; //very small number avoid div-by-zero error  
		double mz1, mz2; 
		int i = 0, j = 0;
		//System.out.println("Matching  " + this.peptide + " with " + s1.peptide);
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("matched " + mz1 + "\t" + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += Math.pow(this.peaks.get(i).getIntensity()
					* s1.peaks.get(j).getIntensity(), 0.5);
				projectedNorm += s1.peaks.get(j).getIntensity();
				i++;
				j++;
			}else{
				j++;
			}
		}
		//System.out.println("this norm: " + this.magnitude());
		//System.out.println("mixture projected norm: " + Math.pow(projectedNorm, 0.5));
		magnitude = magnitude * Math.pow(projectedNorm, 0.5);
		return product/magnitude;
	}
	
	
	public double alpha(Spectrum a, Spectrum b) {		
		double A, B, C, D, E ;		
		C = this.dot(a) ;
		D = a.dot(b) ;
		E = this.dot(b) ;
		A = a.dot(a) ;
		B = b.dot(b) ;
		double alpha = ((B*C)-(D*E) )/ ((A*E)-(C*D));
//		if(alpha > 0){  //there is some strange cases where alpha is < 0
//			return alpha;
//		}else{
//			//System.out.println("warning alpha smaller than zero");
//			return -1*alpha;
//		}
		//we restrict alpha to be smaller than 10, i.e. we require the relative fraction to be no more than 1:10
		if(alpha < 0.1 || alpha > 10){
			//return -1*alpha;
			return 0.1;
		}else{
			if(alpha > 1){
				return 1.0/alpha;
			}else{
				return alpha;
			}
		}
	}
	
	public double maxScore(Spectrum a, Spectrum b, double alpha){
		alpha = Math.pow(alpha, 0.5);  //needed this to take into account of sqrt transformation
		Spectrum answerMix1 = new Spectrum(a, b, alpha, 1) ;
		Spectrum answerMix2 = new Spectrum(b, a, alpha, 1) ;
		//answerMix1.windowFilterPeaks(10, 25);
		//answerMix2.windowFilterPeaks(10, 25);
		Spectrum baseAnswer = new Spectrum(a, b, 1, 0.1); //in case it is a single-match
	    //answerMix1.filterPeaks(100);
	    //answerMix2.filterPeaks(100);
		//answerMix1.sqrtSpectrum();
		//answerMix2.sqrtSpectrum();
	
		//double score1 = (this.dot(answerMix1)) / (this.sumOfPeaks() * answerMix1.sumOfPeaks()) ; 
		//double score2 = (this.dot(answerMix2)) / (this.sumOfPeaks() * answerMix2.sumOfPeaks()) ; 
		double score1 = this.cosineSim(answerMix1);
		double score2 = this.cosineSim(answerMix2);
		double score = this.cosineSim(baseAnswer);
		//double score1 = this.shiftCosineSim(answerMix1);
		//double score2 = this.shiftCosineSim(answerMix2);
	
		if (score1 > score2 && score1 > score)
			return score1 ;
		else if(score2 > score)
			return score2 ;
		else
			return score;
	
	}
	
	
	//shift-version of maxscore, we actually need to try all combination
	/*public double shiftMaxScore(Spectrum a, Spectrum b, double alpha){
	    Spectrum normVa = a.toNormVector();
	    Spectrum normVb = b.toNormVector();
	    
	    Spectrum original = this.toNormVector();
		this.shiftSpectrum(0.5);
		Spectrum rightShift = this.toNormVector();
		this.shiftSpectrum(-1.0); //left shift original spectrum by -0.5 (0 + 0.5 - 1.0) = -0.5
		Spectrum leftShift = this.toNormVector(); 
		this.shiftSpectrum(0.5); //restore the original spectrum
		
		alpha = original.alpha(normVa, normVb);
		double maxOriginal = original.maxScore(normVa, normVb, alpha);
		alpha = rightShift.alpha(normVa, normVb);
		double maxRight = rightShift.maxScore(normVa, normVb, alpha);
		alpha = leftShift.alpha(normVa, normVb);
		double maxLeft = leftShift.maxScore(normVa, normVb, alpha);
		return Math.max(maxLeft,
				Math.max(maxRight, maxOriginal));
	}*/
	
	public double shiftMaxScore(Spectrum a, Spectrum b, double alpha){
	    Spectrum origA = a.toNormVector();
	    a.shiftSpectrum(0.5);
	    Spectrum rightA = a.toNormVector();
	    a.shiftSpectrum(-1.0);
	    Spectrum leftA = a.toNormVector();
	    a.shiftSpectrum(0.5);
	    Spectrum origB = b.toNormVector();
	    b.shiftSpectrum(0.5);
	    Spectrum rightB = b.toNormVector();
	    b.shiftSpectrum(-1.0);
	    Spectrum leftB = b.toNormVector();
	    b.shiftSpectrum(0.5);
	    Spectrum mix = this.toNormVector();
	    
	    double best = 0;
	    alpha = mix.alpha(origA, origB);
		best = Math.max(best, mix.maxScore(origA, origB, alpha));
		alpha = mix.alpha(origA, leftB);
		best = Math.max(best, mix.maxScore(origA, leftB, alpha));
		alpha = mix.alpha(origA, rightB);
		best = Math.max(best, mix.maxScore(origA, rightB, alpha));
		
		alpha = mix.alpha(leftA, origB);
		best = Math.max(best, mix.maxScore(leftA, origB, alpha));
		alpha = mix.alpha(leftA, leftB);
		best = Math.max(best, mix.maxScore(leftA, leftB, alpha));
		alpha = mix.alpha(leftA, rightB);
		best = Math.max(best, mix.maxScore(leftA, rightB, alpha));
		
		alpha = mix.alpha(rightA, origB);
		best = Math.max(best, mix.maxScore(rightA, origB, alpha));
		alpha = mix.alpha(rightA, leftB);
		best = Math.max(best, mix.maxScore(rightA, leftB, alpha));
		alpha = mix.alpha(rightA, rightB);
		best = Math.max(best, mix.maxScore(rightA, rightB, alpha));
		
		return best;
		
	}
	
	//this is an inverse of max score in the sense that rather than dampening the
	//effect of "weaker" spectrum, we actually try to magnify it in the mixture
	//so we pay roughly same attention to both spectrum
	public double inverseMaxScore(Spectrum a, Spectrum b, double alpha){
		//alpha = 10;
		//b.filterPeaks(20);
		Spectrum answerMix = new Spectrum(a, b);
		Spectrum mix = this.toNormVector();
		//mix = this.toNormVector();
		mix.upscale(a, b, alpha);
		//answerMix = answerMix.toNormVector();
		//System.out.println("one part: " + mix.cosineSim(a));
		//System.out.println("second part: " + mix.cosineSim(b));
		//System.out.println("both parts together: " + mix.cosineSim(answerMix));
		//System.out.println(answerMix);
		//System.out.println(mix);
		//System.out.println(a);
		//System.out.println(b);
		return mix.cosineSim(answerMix);
	}
	
	//transforming a mixture spectrum by upscaling one of its component
	//note here we made assumption alpha > 1, need to see cases where it is smaller
	//than one to see if the procedure still hold
	public void upscale(Spectrum s1, Spectrum s2, double alpha){
		//System.out.println(s1);
		//System.out.println(s2);
		double mz1, mz2, mz3; 
		int i = 0, j = 0, k=0;
		double intensity;
		while(i < this.peaks.size() && j < s2.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s2.peaks.get(j).getMass();
			if(k >= s1.peaks.size()){
				mz3 = 3000; //a really big mass
			}else{
				mz3 = s1.peaks.get(k).getMass();
			}
			if(mz1 < mz2){
				i++;
			}else if(mz3 < mz2){
				k++;
			}else if(mz1 == mz2 && mz3 > mz2){
				if(s2.peaks.get(j).getIntensity() > 0.1){
					//System.out.print("s2 intensity: " + s2.peaks.get(j));
					//System.out.print("unscaled intensity: " + this.peaks.get(i));
					this.peaks.get(i).scaleIntensity(alpha+1); //upscale mix intensity
					//System.out.println("scaled intensity: " + this.peaks.get(i).getIntensity());
				}
				i++; 
				j++;
			}else if(mz1 == mz2 && mz1 == mz3){
				if(s2.peaks.get(j).getIntensity() > 0.1){
					intensity = this.peaks.get(i).getIntensity();
					intensity = (intensity - (s1.peaks.get(k).getIntensity()));
					//System.out.print("s2 intensity: " + s2.peaks.get(j));
					//System.out.print("combined intensity: " + this.peaks.get(i));
					//System.out.print("s1 intensity: " + s1.peaks.get(k));
					//System.out.println("subtracted intensity: " + intensity);
					intensity = intensity*(alpha+1);	
					//System.out.println("scaled intensity: " + intensity*(alpha+1));
					if(intensity > 0){
						this.peaks.get(i).setIntensity(this.peaks.get(i).getIntensity() + intensity);
					}
					//System.out.println("final combined intensity: " + this.peaks.get(i));
					
				}
				i++; 
				j++; 
				k++;
			}else{
				j++;
			}
		}
	}
	
	public double dot(Spectrum s1) {
		double product = 0;
		double mz1, mz2; 
		int i = 0, j = 0;
		
		while(i < this.peaks.size() && j < s1.peaks.size()){
			mz1 = this.peaks.get(i).getMass();
			mz2 = s1.peaks.get(j).getMass();
			if(mz1 < mz2){
				i++;
			}else if(mz1 == mz2){
				//System.out.println("Intensity are: " + this.peaks.get(i).getIntensity() + "\t" + s1.peaks.get(j).getIntensity());
				product += Math.pow(this.peaks.get(i).getIntensity(), 2)
					* Math.pow(s1.peaks.get(j).getIntensity(), 2);
				i++;
				j++;
			}else{
				j++;
			}
		}
				
		return product ;
	}
	
	
	//remove peaks that present in both spectrum from this spectrum
	//if magnitue is not the same subtract the intesity of s from this peaks
	//if the resulting peaks is less than zero then the peaks is removed
	public void removeSharePeaks(Spectrum s){
		Iterator<Peak> p1 = this.peaks.iterator();
		Iterator<Peak> p2 = s.getPeak().iterator();
		Peak peak1 = null, peak2 = null;
		if(p1.hasNext() && p2.hasNext()){
			peak1 = p1.next();
			peak2 = p2.next();
		}else{
			return;
		}
		while(p1.hasNext() && p2.hasNext()){
			if(peak1.getMass() < peak2.getMass()){
				peak1 = p1.next();
			}else if(peak1.getMass() == peak2.getMass()){
				p1.remove();
				peak1 = p1.next();
				peak2 = p2.next();
			}else{
				peak2 = p2.next();
			}
		}
		if(peak1.getMass() == peak2.getMass()){ //take care of last element
				p1.remove();
		}
		
	}
	
	public Spectrum removeSharePeaks(Spectrum s, double tolerance){
		int i = 0, j = 0;
		List<Peak> toBeRemoved = new ArrayList<Peak>();
		Spectrum modified = new Spectrum(this);
		while(i < modified.peaks.size() && j < s.peaks.size()){
			Peak m1 = modified.getPeaks().get(i);
			Peak m2 = s.getPeaks().get(j);
			if(m1.getMass() < m2.getMass()){
				i++;
			}
			
			if(m1.getMass() >= m2.getMass()){
				j++;
			}
			
			if(Math.abs(m1.getMass()-m2.getMass()) < tolerance){
				toBeRemoved.add(m1);
			}
		}
		//System.out.println("removing peaks: " + toBeRemoved.size());
		modified.peaks.removeAll(toBeRemoved);
		return modified;
	}
	
	
	public Spectrum project(Spectrum s, double tolerance){
		int i = 0, j = 0;
		List<Peak> toBeAdded = new ArrayList<Peak>();
		Spectrum projected = new Spectrum(this);		
		while(i < projected.peaks.size() && j < s.peaks.size()){
			Peak m1 = projected.getPeaks().get(i);
			Peak m2 = s.getPeaks().get(j);
			if(Math.abs(m1.getMass()-m2.getMass()) < tolerance){
				toBeAdded.add(m1);
				i++;
			}else if(m1.getMass() < m2.getMass()){
				i++;
			}else if(m1.getMass() >= m2.getMass()){
				j++;
			}
			
		}
		//System.out.println("proj peaks: " + toBeAdded.size());
		projected.setPeaks(toBeAdded);
		return projected;
	}
	
	//return the projection as an array, for missing peaks, array has zero value
	//otherwise it is the intensity of the matched peaks
	public double[] projectArray(Spectrum s, double tolerance){
		int i = 0, j = 0;
		Spectrum projected = this;
		Peak dummy = new Peak(0,0);
		Peak match = dummy;
		double[] projection = new double[s.getPeak().size()];
		while(i < projected.peaks.size() && j < s.peaks.size()){
			Peak m1 = projected.getPeaks().get(i);
			Peak m2 = s.getPeaks().get(j);
			//if(Math.abs(m1.getMass()-m2.getMass()) < tolerance){
			if(SWATHUtils.checkMass(m1.getMass(), m2.getMass(), 36, SWATHUtils.PPM)){
				match = m1.getIntensity() > match.getIntensity() ? m1 : match;
				projection[j] += m1.getIntensity();
				i++;
			}else if(m1.getMass() < m2.getMass()){
				i++;
			}else if(m1.getMass() >= m2.getMass()){
				//projection[j] = match.getIntensity();
				match=dummy;
				j++;
			}
			
		}
		//System.out.println(Arrays.toString(projection));
		return projection;
	}
	
	
	
	//removed unfragmented precursors from the spectra
	public void removePrecursors(double tolerance){
		List<Peak> toBeRemoved = new ArrayList<Peak>();
		for(int i = 0; i < this.peaks.size(); i++){
			Peak p = this.peaks.get(i);
			if(Math.abs(p.getMass() - this.parentMass) < tolerance){
				//System.out.println("removing " + p);
				toBeRemoved.add(p);
			}
			
			if(Math.abs(p.getMass() - 
					(this.parentMass - Mass.WATER/this.charge)) < tolerance){
				//System.out.println("removing " + p);
				toBeRemoved.add(p);
			}
			
			if(Math.abs(p.getMass() - 
					(this.parentMass - Mass.NH3/this.charge)) < tolerance){
				//System.out.println("removing " + p);
				toBeRemoved.add(p);
			}
			
			if(Math.abs(p.getMass() - 
					(this.parentMass - 2*Mass.WATER/this.charge)) < tolerance){
				//System.out.println("removing " + p);
				toBeRemoved.add(p);
			}
			
			if(Math.abs(p.getMass() - 
					(this.parentMass - (Mass.WATER+Mass.NH3)/this.charge)) < tolerance){
				//System.out.println("removing " + p);
				toBeRemoved.add(p);
			}
		}
		this.peaks.removeAll(toBeRemoved);
	}
	
	public double residual(Spectrum s){
		Iterator<Peak> p1 = this.peaks.iterator();
		Iterator<Peak> p2 = s.getPeak().iterator();
		Peak peak1 = null, peak2 = null;
		int exp = 4;
//		if(this.isSqrtTrans){
//			exp = 4;
//		}else{
//			exp = 2;
//		}
		double shareIntensity = 0;
		double residIntensity = 0.001; //avoid div-by-zero 
		if(p1.hasNext() && p2.hasNext()){
			peak1 = p1.next();
			peak2 = p2.next();
		}else{
			return shareIntensity/residIntensity;
		}
		double remain;
		//System.out.println("magnitude is: " + this.magnitude());
		while(p1.hasNext() && p2.hasNext()){
			if(peak1.getMass() < peak2.getMass()){
				residIntensity += Math.pow(peak1.getIntensity(),exp);
				//System.out.println("res intensity: " + peak1.getMass() + "\t" + peak1.getIntensity());
				peak1 = p1.next();
			}else if(peak1.getMass() == peak2.getMass()){
				remain = peak1.getIntensity() - peak2.getIntensity();
				if(remain < 0){
					remain = 0;
				}
				shareIntensity += Math.pow(peak1.getIntensity(),exp);
				//System.out.println("share: " + peak2.getMass() + "\t" + peak2.getIntensity());
				//System.out.println("remaining: " + peak1.getMass() + "\t" + remain);
				//residIntensity += Math.pow(remain,exp);
				peak1 = p1.next();
				peak2 = p2.next();
			}else{
				//shareIntensity += Math.pow(peak1.getIntensity(), 2);
				peak2 = p2.next();
			}
		}
		
		if(peak1.getMass() == peak2.getMass()){ //take care of last element
				shareIntensity += Math.pow(peak1.getIntensity(), exp);
				remain = peak1.getIntensity() - peak2.getIntensity();
				if(remain < 0){
					remain = 0;
				}
				//residIntensity += Math.pow(remain, 2);
		}
		if(!p1.hasNext()){
			residIntensity += Math.pow(peak1.getIntensity(), exp);
		}
		
//		if(!p2.hasNext()){
//			shareIntensity += Math.pow(peak2.getIntensity(), 2);
//			//System.out.println("share: " + peak2.getMass() + "\t" + peak2.getIntensity());
//		}
		
		while(p1.hasNext()){
			residIntensity += Math.pow(p1.next().getIntensity(), exp);
		}
		
//		while(p2.hasNext()){
//			shareIntensity += Math.pow(p2.next().getIntensity(), 2);
//		}
			
		//System.out.println("shareIntensity: " + Math.pow(shareIntensity,0.5));
//		System.out.println("s has intensity: " + s.magnitude());
		//System.out.println("resIntensity: " + Math.pow(residIntensity,0.5));
		//double alpha = Math.pow(shareIntensity/residIntensity, 0.5); //initial estimate
		double alpha = Math.pow(residIntensity, 0.5);
		//System.out.println("alpha: " + alpha);
		alpha = alpha*alpha / (1-alpha*alpha);  //reestimate by taking into account that mixture is normalized to one
		alpha =  Math.pow(alpha, 0.5); //thus each component is down-weighted slightly
		if(alpha > 1.0){
			return 1/alpha;
		}else{
			return alpha;
		}
		//return alpha;
	}
	//rather than removing the share peaks, we substract the share
	//peaks from this spectrum, this method is destructive
	public void subtractSharePeaks(Spectrum s){
		Iterator<Peak> p1 = this.peaks.iterator();
		Iterator<Peak> p2 = s.getPeak().iterator();
		Peak peak1 = null, peak2 = null;
		if(p1.hasNext() && p2.hasNext()){
			peak1 = p1.next();
			peak2 = p2.next();
		}else{
			return;
		}
		while(p1.hasNext() && p2.hasNext()){
			if(peak1.getMass() < peak2.getMass()){
				peak1 = p1.next();
			}else if(peak1.getMass() == peak2.getMass()){
				peak1.setIntensity(peak1.getIntensity() - peak2.getIntensity());
				if(peak1.getIntensity() < 0){
					p1.remove();
				}
				peak1 = p1.next();
				peak2 = p2.next();
			}else{
				peak2 = p2.next();
			}
		}
		if(peak1.getMass() == peak2.getMass()){ //take care of last element
				p1.remove();
		}
		
	}
	
	
	public void subtractProjPeaks(Spectrum s, double tolerance){
		Iterator<Peak> p1 = this.peaks.iterator();
		Iterator<Peak> p2 = s.getPeak().iterator();
		Peak peak1 = null, peak2 = null;
		if(p1.hasNext() && p2.hasNext()){
			peak1 = p1.next();
			peak2 = p2.next();
		}else{
			return;
		}
		while(p1.hasNext() && p2.hasNext()){
			if(peak1.getMass() < peak2.getMass()){
				peak1 = p1.next();
			}else if(Math.abs(peak1.getMass() - peak2.getMass()) < tolerance){
				peak1.setIntensity(peak1.getIntensity() - peak2.getIntensity());
			}else{
				peak2 = p2.next();
			}
		}
		if(peak1.getMass() == peak2.getMass()){ //take care of last element
				p1.remove();
		}
		
	}
	
	//filter the spectrum and only keep top n peaks
	public void filterPeaks(int n){
		Vector<Peak> sortedPeakList = new Vector<Peak>();
		sortedPeakList.addAll(this.peaks);
		Collections.sort(sortedPeakList, new peakComparator());
		int i = 0;
		Peak p;
		for(i = 0; i < sortedPeakList.size() - n - 1; i++){
			this.peaks.remove(sortedPeakList.get(i));
		}
	}
	
	public void filterPeaksByMass(double minMass, double maxMass){
		List<Peak> toBeRemoved =  new ArrayList();
		for(int i = 0; i < this.peaks.size(); i++){
			Peak p = this.peaks.get(i);
			if(p.getMass() < minMass){
				toBeRemoved.add(p);
			}else{
				break;
			}

		}
		
		for(int i = this.peaks.size()-1; i > 0; i--){
			Peak p = this.peaks.get(i);
			if(p.getMass() > maxMass){
				toBeRemoved.add(p);
			}else{
				break;
			}
		}
		this.peaks.removeAll(toBeRemoved);
		//System.out.println("removed: " + toBeRemoved);
	}
	
	public void removePeaksInMass(double minMass, double maxMass){
		List<Peak> toBeRemoved =  new ArrayList();
		for(int i = 0; i < this.peaks.size(); i++){
			Peak p = this.peaks.get(i);
			if(p.getMass() > minMass && p.getMass() < maxMass){
				toBeRemoved.add(p);
			}else{
				break;
			}

		}		
		this.peaks.removeAll(toBeRemoved);
		//System.out.println("removed: " + toBeRemoved);
	}
	
	
	//filter spectrum by assuming peak intensity assume normal and keep peaks
	//that are above noise level
	public void filterPeaksByZScore(double  minZ, double lowFract){
		List<Peak> sortedPeakList = new Vector<Peak>();
		List<Peak> lowIntPeaks = new Vector<Peak>();
		List<Peak> noisePeaks = new Vector<Peak>();
		sortedPeakList.addAll(this.peaks);
		Collections.sort(sortedPeakList, new peakComparator());
		int i = 0;
		Peak p;
		int index = (int)Math.ceil(lowFract*sortedPeakList.size());
		for(i = 0; i < index; i++){
			lowIntPeaks.add(sortedPeakList.get(i));
		}
		double mean = this.getMeanIntensity(lowIntPeaks);
		double dev = this.getIntensityDeviation(lowIntPeaks, mean);
		System.out.println("mean: " + mean);
		System.out.println("dev: " + dev);
		for(i = 0; i < sortedPeakList.size(); i++){
			p = this.getPeak().get(i);
			double z = (p.getIntensity() - mean)/dev;
			if(z <= 5){
				noisePeaks.add(p);
			}
		}
		this.peaks.removeAll(noisePeaks);
	}
	
	public void filterPeaksByRankScore(int minR){
		minR = minR - 1; //rank start at -1 by convention;
		List<Peak> sortedPeakList = new Vector<Peak>();
		List<Peak> lowIntPeaks = new Vector<Peak>();
		List<Peak> noisePeaks = new Vector<Peak>();
		for(int i = 0; i < this.peaks.size(); i++){
			Peak p = this.getPeak().get(i);
			if(p.getRank() < minR){
				noisePeaks.add(p);
			}
		}
		int before = this.peaks.size();
		this.peaks.removeAll(noisePeaks);
		System.out.println(this.spectrumName +  "\trank-score filtering\t" + before +"\t" + this.peaks.size());

	}
	
	public void computePeaksZScore(double lowFract){
		List<Peak> sortedPeakList = new Vector<Peak>();
		List<Peak> lowIntPeaks = new Vector<Peak>();
		sortedPeakList.addAll(this.peaks);
		Collections.sort(sortedPeakList, new peakComparator());
		int i = 0;
		Peak p;
		int index = (int)Math.ceil(lowFract*sortedPeakList.size());
		//int index = sortedPeakList.size() - 500;
		//index = index < 2 ?  (int)Math.ceil(lowFract*sortedPeakList.size()) : index;
		for(i = 0; i < index; i++){
			lowIntPeaks.add(sortedPeakList.get(i));
		}
		if(lowIntPeaks.size() < 2){
			return;
		}
		double mean = this.getMeanIntensity(lowIntPeaks);
		double median = lowIntPeaks.get((int)(lowIntPeaks.size()/2.0)+1).getIntensity();
		double dev = this.getIntensityDeviation(lowIntPeaks, mean);
		System.out.println("Scan: " + this.scanNumber + "\tmean:\t" + mean +  "\tmedian\t" + median + "\tdev:\t" + dev);
		for(i = 0; i < sortedPeakList.size(); i++){
			p = this.getPeak().get(i);
			//double z = (p.getIntensity() - mean)/dev;
			double z = p.getIntensity()/median*100;
			p.setRank((int)z);
		}
	}
	
	public void computeAnnotatedZScore(Spectrum s, double tolerance){
		Peptide pep = new Peptide(s.peptide, s.charge);
		String[] prefixIons = {"b", "b-H20", "b-NH3", "b(iso)"};//"a", "a-H20", "a-NH3"};
		String[] suffixIons = {"y", "y-H20", "y-NH3", "y(iso)"};
		TheoreticalSpectrum.prefixIons = prefixIons;
		TheoreticalSpectrum.suffixIons = suffixIons;
		TheoreticalSpectrum t = new TheoreticalSpectrum(pep);
		//TheoreticalSpectrum.addIsotopicPeaks(t, 1);
		SimpleMatchingGraph g = t.getMatchGraph(s, tolerance);
		List<Peak> annotated = new ArrayList();
		List<Peak> unAnnotated = new ArrayList();
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				unAnnotated.add(p);
			}else{
				annotated.add(p);
			}
		}
		Collections.sort(unAnnotated, PeakIntensityComparator.comparator);
		double mean = 0.0;
		double topFract = 0.01;
		double minZ = 150;
		int length = unAnnotated.size();
		int count=0;
		for(int i = length-1; i > (int)(length*topFract); i--){
			mean+=unAnnotated.get(i).getIntensity();
			count++;
		}
		
		mean = mean / count;
		
		s.mergePeaks(s, 0.05);
		g = t.getMatchGraph(s, tolerance);
		unAnnotated.clear();
		annotated.clear();
		for(Iterator<LabelledPeak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = g.getNeighbors(p);
			if(neighs.size() == 0){
				unAnnotated.add(p);
			}else{
				annotated.add(p);
			}
		}
		
		int bigSignal = 0;
		for(int i = 0; i < annotated.size(); i++){
			double z = (annotated.get(i).getIntensity() / mean)*100;
			annotated.get(i).setRank((int)z);
			if(z > minZ){
				bigSignal++;
				//System.out.println("big signal: " + annotated.get(i));
			}
		}
		
		int bigNoise = 0;
		for(int i = 0; i < unAnnotated.size(); i++){
			double z = (unAnnotated.get(i).getIntensity() / mean)*100;
			unAnnotated.get(i).setRank((int)z);
			if(z > minZ){
				bigNoise++;
				//System.out.println("big nosie: " + unAnnotated.get(i));
			}
		}
		
		System.out.println(s.spectrumName + "\t" + s.peptide + "\t" + s.charge + "\ttotal peaks:\t" 
				+ annotated.size() + "\t" + unAnnotated.size() + "\tmean-noise level:\t" 
				+ mean + "\tabove-" +minZ + "%:\t" + bigSignal + "\t" + bigNoise);
	}
	
	
	
	public void filterPeaksByIntensity(double minInt){
		List<Peak> noiseList = new ArrayList<Peak>();
		for(int i = 0; i < this.peaks.size(); i++){
			Peak p = this.peaks.get(i);
			if(p.getIntensity() < minInt){
				noiseList.add(p);
			}
		}
		int before = this.getPeaks().size();
		this.peaks.removeAll(noiseList);
		//System.out.println(this.spectrumName + "\tintensity-filtering\t" +before + "\t" +  this.peaks.size());
	}
	
	public void filterFlatPeaks(double minFraction){
		Map<Double, List<Peak>> intensityFreq = new HashMap<Double, List<Peak>>();
		List<Peak> pList = null;
		List<Peak> noiseList = new ArrayList<Peak>();
		for(int i = 0; i < this.getPeak().size(); i++){
			Peak p = this.getPeak().get(i);
			if(intensityFreq.containsKey(p.getIntensity())){
				pList = intensityFreq.get(p.getIntensity());
			}else{
				pList = new ArrayList<Peak>();
			}
			pList.add(p);
			intensityFreq.put(p.getIntensity(), pList);
		}
		//this.peaks.removeAll(noisePeaks);
		for(Iterator<Double> it = intensityFreq.keySet().iterator(); it.hasNext();){
			double key = it.next();
			int count = intensityFreq.get(key).size();
			if((double)count/this.peaks.size() > minFraction){
				noiseList.addAll(intensityFreq.get(key));
			}
		}
		this.peaks.removeAll(noiseList);
	}
	
	//for spectrum after sqrt transformation
	public int explainedIntensity(double percent){
		Vector<Peak> sortedPeakList = new Vector<Peak>();
		sortedPeakList.addAll(this.peaks);
		Collections.sort(sortedPeakList, new peakComparator());
		double mag = this.magnitude();
		mag = mag*mag;
		int i = 0;
		Peak p;
		double current = 0.0;
		for(i = sortedPeakList.size()-1; i > 0;  i--){
			current += sortedPeakList.get(i).getIntensity()
				*sortedPeakList.get(i).getIntensity();
			if(current / mag > percent){
				return sortedPeakList.size() - i;
			}
		}
		return sortedPeakList.size();
	}
	
	//for spectrum before sqrt
	public int explainedIntensity2(double percent){
		Vector<Peak> sortedPeakList = new Vector<Peak>();
		sortedPeakList.addAll(this.peaks);
		Collections.sort(sortedPeakList, new peakComparator());
		double mag = this.sumMagnitude();
		int i = 0;
		Peak p;
		double current = 0.0;
		for(i = sortedPeakList.size()-1; i > 0;  i--){
			current += sortedPeakList.get(i).getIntensity();
			if(current / mag > percent){
				return sortedPeakList.size() - i;
			}
		}
		return sortedPeakList.size();
	}
	
	public int intensePeakCount(double threshold){
		int count = 0;
		for(int i = 0; i < this.peaks.size(); i++){
			if(peaks.get(i).getIntensity() > threshold){
				count++;
			}
		}
		return count;
	}
	//window based filtering, for every peak m, we look in
	//the range +/- deltaM from m and see if peak m is ranked
	//top n peaks in the window if it is, the peak is kept otherwise it isremoved
	public void windowFilterPeaks(int n, double deltaM){
		int i = 0;
		Peak p;
		Vector<Peak> toBeRemove = new Vector();
		Vector<Peak> neighbors = new Vector();
		for(int j = 0; j < this.peaks.size(); j++){
			//System.out.println("processing peak " + j);
			p = this.peaks.get(j);
			//System.out.println("we had " + neighbors.size() + " neighbors");
			removeNonNeighborPeaks(deltaM, neighbors, p);
			//System.out.println("we have " + neighbors.size() + " remaining neighbors");
			i = addNeighborPeaks(deltaM, neighbors, p, i);
			//System.out.println("we end up with " + neighbors.size() + " neighbors this time");
			Vector sortedNeighbors = new Vector(neighbors);
			Collections.sort(sortedNeighbors, new peakComparator());
			//System.out.println("rank is: " + getPeakRank(sortedNeighbors, p));
			if(getPeakRank(sortedNeighbors, p) > n){
				toBeRemove.add(p);
			}
		}
		//System.out.println("We start with " + peaks.size() + " peaks");
		this.peaks.removeAll(toBeRemove);
		//System.out.println("After window-filtering we have: " + peaks.size() + " peaks");		
	}
	
	public void windowFilterPeaks2(int n, double deltaM){
		int current = 0, left = 0, right = 0;
		Collection<Peak> toBeRemove = new ArrayList();
		TreeSet<Peak> neighs = new TreeSet(PeakIntensityComparator.comparator);
		while(current < this.peaks.size()){
			Peak p = this.peaks.get(current);
			for(int i = left; i < current; i++){
				Peak smaller = this.peaks.get(i);
				if(p.getMass() - smaller.getMass() > deltaM){
					neighs.remove(smaller);
					left = i;
				}
			}
			for(int i = right+1;  i < this.peaks.size(); i++){
				Peak bigger = this.peaks.get(i);
				bigger.setIntensity(bigger.getIntensity()+Math.random()*0.000001);
				if(bigger.getMass() - p.getMass()  <= deltaM){
					neighs.add(bigger);
				}else{
					right = i-1;
					break;
				}
			}
			//System.out.println("current list size: " + neighs.size());
			Iterator<Peak> it = neighs.descendingIterator();
			for(int i = 0; it.hasNext(); i++){
//				//System.out.println("i is: " +i);
				Peak p2 = it.next();
				//System.out.println(p2);
//				//System.out.println("peak " + p2);
				if(p == p2){
					break;
				}
				if(i > n){
					toBeRemove.add(p);
					break;
//					//System.out.println("remove " + p);
				}
			}
			//System.out.println();
			current++;
		}
		int before = this.getPeak().size();
		this.peaks.removeAll(toBeRemove);
		//System.out.println(this.spectrumName + "\twindow-filtering:\t" + before + "\t" + this.getPeak().size());
	}
	
	public void windowFilterAndRank(int n, double deltaM, int topN){
		int i = 0;
		Peak p;
		Vector<Peak> toBeRemove = new Vector();
		Vector<Peak> neighbors = new Vector();
		for(int j = 0; j < this.peaks.size(); j++){
			p = this.peaks.get(j);
			//System.out.println("we had " + neighbors.size() + " neighbors");
			removeNonNeighborPeaks(deltaM, neighbors, p);
			//System.out.println("we have " + neighbors.size() + " remaining neighbors");
			i = addNeighborPeaks(deltaM, neighbors, p, i);
			//System.out.println("we end up with " + neighbors.size() + " neighbors this time");
			Vector sortedNeighbors = new Vector(neighbors);
			Collections.sort(sortedNeighbors, new peakComparator());
			//System.out.println("rank is: " + getPeakRank(sortedNeighbors, p));
			if(getPeakRank(sortedNeighbors, p) > n){
				toBeRemove.add(p);
			}
		}
		this.peaks.removeAll(toBeRemove);
		Collections.sort(toBeRemove, PeakIntensityComparator.comparator);
		this.computePeakRank();
		int count=1;
		int toBeInsert = topN-this.peaks.size()+1;
		for(int k = toBeRemove.size()-1; k > toBeRemove.size() - toBeInsert && k > 0; k--){
			toBeRemove.get(k).setRank(this.peaks.size()+count);
			this.peaks.add(toBeRemove.get(k));
			count++;
		}
		Collections.sort(this.peaks, PeakMassComparator.comparator);
		
	}
	
	public List<Peak> topWindowPeaks(int n, double deltaM){
		int i = 0;
		Peak p;
		List<Peak> toBeKept = new Vector<Peak>();
		List<Peak> neighbors = new Vector<Peak>();
		for(int j = 0; j < this.peaks.size(); j++){
			p = this.peaks.get(j);
			//System.out.println("we had " + neighbors.size() + " neighbors");
			removeNonNeighborPeaks(deltaM, neighbors, p);
			//System.out.println("we have " + neighbors.size() + " remaining neighbors");
			i = addNeighborPeaks(deltaM, neighbors, p, i);
			//System.out.println("we end up with " + neighbors.size() + " neighbors this time");
			List<Peak> sortedNeighbors = new Vector<Peak>(neighbors);
			Collections.sort(sortedNeighbors, new peakComparator());
			//System.out.println("rank is: " + getPeakRank(sortedNeighbors, p));
			if(getPeakRank(sortedNeighbors, p) <= n){
				toBeKept.add(p);
			}
		}
		System.out.println("We start with " + peaks.size() + " peaks");
		System.out.println("After window-filtering we have: " + toBeKept.size() + " peaks");	
		return toBeKept;
	}
	
	public Map<Peak, Integer> getPeakRank(){
		List<Peak> sortedList = new Vector<Peak>();
		sortedList.addAll(this.peaks);
		Collections.sort(sortedList, new peakComparator());
		Map<Peak, Integer> peakRankMap = new HashMap<Peak, Integer>();
		Peak p;
		for(int i = 0, size = sortedList.size(); i < size; i++){
			p = sortedList.get(i);
			peakRankMap.put(p, new Integer(i));
		}
		return peakRankMap;
	}
	
	public void computePeakRank(){
		List<Peak> sortedList = new Vector<Peak>();
		sortedList.addAll(this.peaks);
		Collections.sort(sortedList, new peakComparator());
		Map<Peak, Peak> peakMap = new HashMap<Peak, Peak>();
		Peak p;
		
		for(int i = 0, size = sortedList.size(); i < size; i++){
			p = sortedList.get(i);
			p.setRank(size - i);
		}
	}
	
	public void computePeakSector(){
		int count1=0, count2=0;
		for(int i = 0, size = this.peaks.size(); i < size; i++){
			Peak p = this.peaks.get(i);
			if(p.getMass() / (this.parentMass*this.charge) < 0.5){
				count1++;
				p.setSector(0);
			}else{
				count2++;
				p.setSector(1);
			}
		}
		System.out.println("sector 0 total: " + count1 + " sector 1 total: " + count2);
	}
	
	public List<Peak> getSortedPeaks(){
		List<Peak> sortedList = new Vector<Peak>();
		sortedList.addAll(this.peaks);
		Collections.sort(sortedList, new peakComparator());
		return sortedList;
	}
	
	public List<Peak> getTopPeaks(int N){
		List<Peak> sortedList = getSortedPeaks();
		int begin = sortedList.size() - N;
		begin = begin < 0 ? 0 : begin;
		return sortedList.subList(begin, sortedList.size());
	}
	
	//require peaks at least certain mass difference aprt
	public List<Peak> getTopPeaks(int N, double minMassDiff){
		List<Peak> sortedList = getSortedPeaks();
		List<Peak> ret = new ArrayList();
		int addCount = 0;
		int index = sortedList.size()-1;
		while(index >= 0 && addCount < N){
			Peak toBeAdd = sortedList.get(index);
			boolean canAdd = true;
			for(int i = 0; i < ret.size(); i++){
				if(Math.abs(ret.get(i).getMass() - toBeAdd.getMass()) < minMassDiff){
					canAdd = false;
					break;
				}
			}
			if(canAdd){
				ret.add(toBeAdd);
				addCount++;
			}
			index--;
		}
		return ret;
	}
	
	public List<Peak> getTopPeaks(int N, List<Peak> exclusionList){
		List<Peak> sortedList = getSortedPeaks();
		List<Peak> out = new ArrayList<Peak>();
		for(int i = sortedList.size() - N; i > 0 && out.size() < N; i--){
			Peak p = sortedList.get(i);
			if(!exclusionList.contains(p)){
				out.add(p);
			}
		}
		return out;
	}
	
	public List<Peak> getProjectedTopPeaks(TheoreticalSpectrum t, int N){
		List<Peak> sortedList = getSortedPeaks();
		SimpleMatchingGraph g = t.getMatchGraph(this, 0.5);
		Set<Peak> unMatchedPeaks = new HashSet<Peak>();
		Set<Peak> actuals = g.vertexSet(1);
		List<Peak> topNPeak = new ArrayList();
		for(Iterator<Peak> it = actuals.iterator(); it.hasNext();){
			Peak curr = it.next();
			List neighbors = g.getNeighbors(curr);
			if(neighbors.size() == 0){
				unMatchedPeaks.add(curr);
			}
		}
		for(int i = sortedList.size()-1, count = 0; i > 0 && count < N; i--){
			Peak curr = sortedList.get(i);
			if(unMatchedPeaks.contains(curr)){
				topNPeak.add(curr);
				count++;
			}
		}
		return topNPeak;
	}
	
	//assume peaks are sorted by mass
	private int addNeighborPeaks(double deltaM, Collection<Peak> v, Peak p, int beginInd){
		boolean beginFlag = false;
		for(int i = beginInd; i < this.peaks.size(); i++){
			if(Math.abs(peaks.get(i).getMass() - p.getMass()) <= deltaM){
				v.add(peaks.get(i));
				beginFlag = true;
			}else if(beginFlag){
				return i;
			}else if(beginFlag && i == this.peaks.size()-1){
				return i+1;
			}
		}
		return beginInd;
		
	}
	
	private void removeNonNeighborPeaks(double deltaM, Collection<Peak> v, Peak p){
		Iterator<Peak> it = v.iterator();
		while(it.hasNext()){
			if(Math.abs(it.next().getMass() - p.getMass()) > deltaM){
				it.remove();
			}
		}
	}
	
	private int getPeakRank(List<Peak> l, Peak p){
		for(int i = 0; i < l.size(); i++){
			if(l.get(i).equals(p)){
				return l.size()-i;
			}
		}
		return -1;
	}
	
	//return the minimum number of peaks that exaplained
	//a percentage of total intensity, this way we do not
	//account for very low-intensity peaks that are likely
	//to be noise rather than true signal
	public int numberOfPeaks(double percent){
		Vector<Peak> sortedPeakList = new Vector<Peak>();
		sortedPeakList.addAll(this.peaks);
		Collections.sort(sortedPeakList, new peakComparator());
		//double total = this.magnitude();
		double total = 0;
		int i = 0, count = 0;
		for(i = sortedPeakList.size()-1; i >= 0; i--){
			total += sortedPeakList.get(i).getIntensity();
		}
		double current = 0;
		for(i = sortedPeakList.size()-1; i >= 0; i--){
			current += sortedPeakList.get(i).getIntensity();
			 //*sortedPeakList.get(i).getIntensity();
			if(current / total >= percent){
				//System.out.println("percent is: " + m/total);
				return ++count;
			}else{
				count++;
			}
			
		}
		return count;
	}
	
	//printout the spectrum for us to see, using default MGF format
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("BEGIN IONS\n");
		//if(peptide != null ){
		//	sb.append("TITLE=" + this.peptide + "\n");
		//}else{
			sb.append("TITLE=" + this.spectrumName + "\n");
			//sb.append("TITLE=NIST spectral library entry: " + this.peptide.split("\\.")[0]  + "\n");
		//}
		if(charge >= 0){	
			sb.append("CHARGE=" + charge + "+\n");
		}else{
			sb.append("CHARGE=" + charge + "-\n");
		}
		sb.append("PEPMASS=" + parentMass + "\n");
		if(this.peptide.equals("DUMMYSEQ")){
			//no annotatin to be print
		}else if(this.peptide != null && this.peptide.contains(" & ")){
			//String[] peps = this.peptide.split(" & ");
			//sb.append("PEPSEQ=" + peps[0].split("\\.")[0] + "\n");
			//sb.append("PEPSEQ=" + peps[1].split("\\.")[0] + "\n");
			sb.append("SEQ="+this.peptide+"\n");
		}else{
			sb.append("SEQ=" + peptide + "\n");
			//sb.append("SEQ=" + peptide.split("\\.")[0] + "\n");
		}
		//sb.append("PEPEXP=" + 0.00000001+ "\n"); //generic expect value
		//sb.append("PEPMOD=" + modMass + "@" + modPos + "\n");
		//sb.append("GPMp\n");	
		if(this.scanNumber > 0){
			sb.append("SCAN="+this.scanNumber+"\n");
		}
		for(int i = 0; i < peaks.size(); i++){
			sb.append(peaks.get(i)+"\n");
		}
		sb.append("END IONS\n");
		return sb.toString();
	}
	
	public static void testZScoreFilter(){
		String spectrumFile = "../mixture_linked/msdata/gringar/swath_development/14304_EIF4A2_SWATH-1x.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		Spectrum s = reader.getSpectrum(21785);
		System.out.println("Number of peaks before: " + s.getPeak().size());
		//s.filterPeaksByZScore(3, 0.9);
		s.filterFlatPeaks(0.1);
		System.out.println("Number of peaks after: " + s.getPeak().size());
		System.out.println(s);
	}
	
	public static void main(String[] args){
		//testReadMGF();
		//testtoVector();
		//testVectorRep();
		//testMixSpect();
		//testRemoveSharePeak();
		//testPeakCount();
		//testfilterPeak();
		//testWindowFilterPeak();
		//testShiftCosine();
		testZScoreFilter();
	}
	//the following are simple test cases for various functions
	public static void testReadMGF(){
		String filename = "testspectrum.mgf";
		Spectrum msms = new Spectrum();
		msms.readSpectrumFromMGF(filename);
		System.out.print(msms);
	}
	
	public static void testtoVector(){
		System.out.println("generating vectored spectrum");
		String filename = "testspectrum.mgf";
		Spectrum msms = new Spectrum();
		msms.readSpectrumFromMGF(filename);
		Spectrum s = msms.toVector(50, 100, 500);
		System.out.println(s);
	}
	
	public static void testVectorRep(){
		String filename = "testspectrum.mgf";
		Spectrum msms = new Spectrum();
		msms.readSpectrumFromMGF(filename);
		Spectrum msms2 = new Spectrum();
		msms2.readSpectrumFromMGF(filename);
		System.out.println("Comparing self to self: "
				+ msms2.cosineSim(msms));
	
		msms = msms.toNormVector(0.5, 100, 500);
		msms2 = msms2.toNormVector(0.5, 100, 500);
		System.out.println("Comparing self to self: "
				+ msms2.cosineSim(msms));
		System.out.println("scaling one copy of self");
		msms.scaleSpectrum(0.03);
		System.out.println("Comparing self to self: "
				+ msms2.cosineSim(msms));
		System.out.println();
	
	
	}
	
	public static void testMixSpect(){
		System.out.println("generating mix spectrum");
		String filename = "testspectrum.mgf";
		Spectrum s1 = new Spectrum();
		s1.readSpectrumFromMGF("testspectrum.mgf");
		Spectrum s2 = new Spectrum();
		s2.readSpectrumFromMGF("testspectrum2.mgf");
		Spectrum msms12 = new Spectrum(s1, s2);
		s1=s1.toNormVector();
		s2=s2.toNormVector();
		System.out.println(s1);
		System.out.println(s2);
		System.out.println(msms12);
		System.out.println("cosine: " + s2.cosineSim(msms12));
		System.out.println("projected: " + s2.projectedCosine(msms12));
		System.out.println("cosine: " + s1.cosineSim(msms12));
		System.out.println("projected: " + s1.projectedCosine(msms12));
		System.out.println("distinct: " + s1.cosineSim(s2));
		System.out.println("distinct projected: " + s1.projectedCosine(s2));
	}
	
	public static void testRemoveSharePeak(){
		String filename = "testspectrum.mgf";
		Spectrum msms = new Spectrum();
		msms.readSpectrumFromMGF(filename);
		Spectrum msms2 = new Spectrum();
		msms2.readSpectrumFromMGF(filename);
		System.out.println("original spectrum");
		System.out.println(msms);
		System.out.println("removing self against self");
		msms.removeSharePeaks(msms2);
		System.out.println(msms);
		msms.readSpectrumFromMGF(filename);
		msms.peaks.get(0).setMoz(20);  //manually changing the first peak 
		System.out.println(msms);
		System.out.println("removing all but the first peak");
		msms.removeSharePeaks(msms2);
		System.out.println(msms);
	}
	
	public static void testPeakCount(){
		String filename = "testspectrum.mgf";
		Spectrum msms = new Spectrum();
		msms.readSpectrumFromMGF(filename);
		System.out.println(msms);
		System.out.println("total magnitude: " + msms.magnitude());
		for(double i = 0; i < 11; i++){
			System.out.println("number of peaks at" +  (i*10) +  "%: " + msms.numberOfPeaks(i/10));
		}
		
	}
	
	public static void testfilterPeak(){
		String filename = "testspectrum.mgf";
		Spectrum msms = new Spectrum();
		msms.readSpectrumFromMGF(filename);
		System.out.println(msms);
		System.out.println("total magnitude: " + msms.magnitude());
		for(int i = 20; i > 0; i--){
			msms.filterPeaks(i);
			System.out.println(msms);
		}
	}
	
	public static void testWindowFilterPeak(){
		String filename = "testspectrum.mgf";
		Spectrum msms = new Spectrum();
		msms.readSpectrumFromMGF(filename);
		System.out.println("before window filter: ");
		System.out.println(msms);
		msms.windowFilterPeaks(2, 50);
		System.out.println("after window filter: ");
		System.out.println(msms);
	}
	
	public static void testShiftCosine(){
		String filename = ".\\mixture_compressed\\new80min.mgf";
		String fileMix = ".\\mixture_compressed\\new2min.mgf";
        SpectrumLib lib1 = new SpectrumLib(filename, "MGF");
		SpectrumLib mixlib = new SpectrumLib(fileMix, "MGF");
		Spectrum s = lib1.getSpectra("spec_3182.dta..1").get(0);
		Spectrum m = mixlib.getSpectra("spec_186.dta..1").get(0);
		Spectrum sCopy = s;
		System.out.println("shift cosine: " + s.shiftCosineSim(m));
		System.out.println("shift cosine: " + s.shiftCosineSim(m));
		System.out.println("shift cosine: " + m.shiftCosineSim(s));
		System.out.println("");
		System.out.println("shift cosine with self: " + s.shiftCosineSim(sCopy));
		System.out.println("shift cosine with self: " + sCopy.shiftCosineSim(s));
		System.out.println("");
		s = s.toNormVector();
		sCopy = sCopy.toNormVector();
		System.out.println("cosine: " + s.cosineSim(sCopy));
		System.out.println("cosine: " + sCopy.cosineSim(s));

	}
	//a simple comparator that compare the intensity of peaks
	private class peakComparator implements Comparator<Peak>{ 
		public int compare(Peak p0, Peak p1) {
			if(p0.getIntensity()> p1.getIntensity()){
				return 1;
			}else if(p0.getIntensity() == p1.getIntensity()){
				return 0;
			}else{
				return -1;
			}
		}
	}
}


 