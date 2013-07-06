package SeqDB;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.Spectrums.Mass;
import org.Spectrums.Peptide;
import org.Spectrums.PeptideLite;

import sequences.FastaSequence;

public class PeptideUtils {
	//reduce redundant peptides in database
	public static List<Peptide> generatePeptide(List<PeptideLite> pList, FastaSequence seq, double tolerance){
		List<Peptide> peptides = new ArrayList();
		Set<String> strs = new HashSet();
		for(int i = 0; i < pList.size(); i++){
			PeptideLite pep = pList.get(i);
			String s = seq.getSubsequence(pep.getBeginInd(), pep.getEndInd()+1);
			//System.out.println("string is : " + s);
			if(!strs.contains(s)){
				Peptide p = new Peptide(s+"."+1);
				peptides.add(p);
				strs.add(s);
			}
		}
		return peptides;
	}
	
	public static List<Peptide> generatePeptide(List<PeptideLite> pList, FastaSequence seq, double fromMass, double toMass){
		return generatePeptide(pList, seq, fromMass, toMass, 1);
	}
	public static List<Peptide> generatePeptide(List<PeptideLite> pList, FastaSequence seq, double fromMass, double toMass, int C13){
		List<Peptide> peptides = new ArrayList();
		Set<String> strs = new HashSet();
		fromMass = fromMass + Mass.WATER + Mass.PROTON_MASS;
		toMass = toMass + Mass.WATER + Mass.PROTON_MASS;
		double offset = Mass.C13 - Mass.C12;
		for(int i = 0; i < pList.size(); i++){
			PeptideLite pep = pList.get(i);
			String s = seq.getSubsequence(pep.getBeginInd(), pep.getEndInd()+1);
			//System.out.println("string is : " + s);
			if(!strs.contains(s)){
				Peptide p = new Peptide(s+"."+1);
				double parentMass = p.getParentmass();
				for(int c = 0; c <= C13; c++){
					double fromMassIso = fromMass - c*offset;
					double toMassIso = toMass - c*offset;
					if(parentMass > fromMassIso && parentMass < toMassIso){
						peptides.add(p);
						strs.add(s);
					}
				}
			}
		}
		return peptides;
	}
	
	//compute mass of a peptide from spectrum info
	public static double getPeptideMass(double precursor, int charge){
		return precursor*charge - Mass.PROTON_MASS*charge - Mass.WATER; 
	}

}
