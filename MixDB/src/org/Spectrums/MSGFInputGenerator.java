package org.Spectrums;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

//generate MSGF output for spectral library searches



public class MSGFInputGenerator {
	public static void generateInputs(){
		String libraryPepFile = "..\\mixture_linked\\spectral_library\\human_speclib_peptide_plus_spectrast_decoy_withCorrectedPTM.txt";
		String spectrumFile = "..\\mixture_linked\\human_heck_data\\data\\090121_NM_Trypsin_24.mzXML";
		List<String> peptides = Utils.FileIOUtils.createListFromFile(libraryPepFile);
		List<Peptide> pepList = new ArrayList<Peptide>();
		int massOffCount = 0, decoyOffCount=0;
		for(int i = 0; i < peptides.size(); i++){
			//System.out.println("line is: " + peptides.get(i));
			String[] tokens = peptides.get(i).split("\\s+");
			String[] tokens2 = tokens[1].split("[\\/]");
			Peptide p = new Peptide(tokens2[0].substring(2, tokens2[0].length()-2), Integer.parseInt(tokens2[1]));
			//System.out.println("peptide is: " + p + "\t" + p.getCharge() + "\t" + p.getParentmass());
			if(Math.abs(Double.parseDouble(tokens[2]) - p.getParentmass()) > 0.01){
				if(tokens[3].equals("decoy")){
					decoyOffCount++;
				}
				massOffCount++;
			}
			p.setParentmass(Double.parseDouble(tokens[2]));
			//System.out.println("line is: " + peptides.get(i));
			//System.out.println("peptide is " + p +"\t" + p.getParentmass() + "\t" + p.getCharge());
			pepList.add(p);
			
			if(tokens[3].equals("decoy")){
				p.setDecoy(true);
			}
		}
		//System.out.println("Finish processing library, found mass error count: " + massOffCount + "\tdecoy mass of count: " + decoyOffCount);
		CandidateSpectrumLibFactory factory = new CandidateSpectrumLibFactory();
		factory.indexPeptideByParentMass(0.5, pepList);
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		String[] tokens = spectrumFile.split("\\\\");
		System.out.println("#SpectrumFile\tTitle\tScan#\tAnnotation\tCharge\tProtein");
		for(;reader.hasNext();){
			Spectrum s = reader.next();
			List<Peptide> candidates = factory.getCandidateByMass(s.parentMass, 3.0);
			for(Iterator<Peptide> it = candidates.iterator(); it.hasNext();){
				Peptide candidate = it.next();
				String protein = "Target_protein";
				if(candidate.isDecoy()){
					protein = "Decoy_protein";
				}
				System.out.println(tokens[tokens.length-1]+"\t"+"title" + "\t" + s.scanNumber + "\t*." + candidate + ".*\t" + candidate.getCharge() + "\t" + protein);
			}
		}
	}
	
	public static void generateInputFromMGFLibrary(){
		String libraryFile = "..\\mixture_linked\\spectral_library\\NIST_yeast_IT_v3.0_2009-05-04_7AA.mgf";
		LargeSpectrumLibIterator iter = new LargeSpectrumLibIterator(libraryFile);
		int counter = 0;
		System.out.println("#SpectrumFile\tTitle\tScan#\tAnnotation\tCharge\tProtein");
		while(iter.hasNext()){
			Spectrum current = (Spectrum)iter.next();
			Peptide p = new Peptide(current.peptide);
			System.out.println(getFileName(libraryFile)+"\t"+"title" + "\t" + counter + "\t*." + p + ".*\t" + current.charge + "\t" + "Protein");
			counter++;
		}
		
	}
	
	public static String getFileName(String pathName){
		String[] tokens = pathName.split("\\\\|/");
		return tokens[tokens.length-1];
	}
	
	public static void processPTM(Peptide p){
		System.out.println("#SpectrumFile\tTitle\tScan#\tAnnotation\tCharge\tProtein");
		for(int i = 0; i < p.getPos().length; i++){
			char residue = p.getPeptide().charAt(p.getPos()[i]-1);
			double aaMass = Mass.getAAMass(residue);
			p.getPtmmasses()[i] = p.getPtmmasses()[i]-aaMass;
		}
	}
	
	public static void main(String[] args){
		generateInputs();
		//generateInputFromMGFLibrary();
	}
}
