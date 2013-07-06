package mixdb;

import java.util.Iterator;
import java.util.Map;

import org.Spectrums.LinkedPeptidePeakScoreLearner;
import org.Spectrums.MixturePeak;
import org.Spectrums.MixturePeakScoreLearner;
import org.Spectrums.Peak;
import org.Spectrums.Peptide;
import org.Spectrums.Spectrum;

public class LinkedScorerAdapter {
	private Map<String, LinkedIonType> typeMap;
	private IonTypeMapper ionTypeMap;
	private LinkedPeptidePeakScoreLearner peakComp;
	private double[][] table; //the new scoring table
	private double[][] errors;
	private boolean DEBUG=false;
	
	public LinkedScorerAdapter(Map<String, LinkedIonType> typeMap, IonTypeMapper ionMap, LinkedPeptidePeakScoreLearner pComp){
		this.typeMap = typeMap;
		this.ionTypeMap = ionMap;
		this.peakComp = pComp;
		getScoringTable();
		getErrorTable();
	}
	
	private void getScoringTable(){
		int[] ranks = this.peakComp.getRankInterval();
		this.table = new double[typeMap.values().size()][ranks.length];
		Peptide[] dummyPeps = {new Peptide("AAAAAAAAA.1")};
		int counter = 0;
		for(Iterator<LinkedIonType> it = typeMap.values().iterator(); it.hasNext();){
			LinkedIonType type = (LinkedIonType)it.next();
			Spectrum s = new Spectrum();
			s.charge = type.getMixtureCharge();
			MixturePeptideType pType = (MixturePeptideType)type.getPType(); 
			Peptide p = dummyPeps[pType.getLength()];
			p.setCharge((short)pType.getPepCharge());
			p.setLinkedPos(3);
			MixturePeak dummyPeak = null;
			if(type.getDirection() == SimpleIonType.PREFIX){
				if(type.isLinked()){	
					dummyPeak = new MixturePeak(0, 0, type.getType(), (short)(p.getLinkedPos()+1), (short)type.getCharge(), pType.getPeptideAbundance()); //make it linked
				}else{
					dummyPeak = new MixturePeak(0, 0, type.getType(), (short)(p.getLinkedPos()-1), (short)type.getCharge(), pType.getPeptideAbundance());  //make it unlinked
				}
			}
			if(type.getDirection() == SimpleIonType.SUFFIX){
				if(type.isLinked()){	
					dummyPeak = new MixturePeak(0, 0, type.getType(), (short)(p.getPeptide().length()-p.getLinkedPos()+1), (short)type.getCharge(), pType.getPeptideAbundance()); //make it linked
				}else{
					dummyPeak = new MixturePeak(0, 0, type.getType(), (short)(p.getPeptide().length()-p.getLinkedPos()), (short)type.getCharge(), pType.getPeptideAbundance());  //make it unlinked
				}
			}
			dummyPeak.setParent(s);
			dummyPeak.setPep(p);
			Peak dummyPeak2 = new Peak(0,0);
			for(int i = 0; i < ranks.length-1; i++){
				//System.out.println("rank is: " + ranks[i]);
				dummyPeak2.setRank(ranks[i]);
				double score = peakComp.getScore(dummyPeak, dummyPeak2);
				table[this.ionTypeMap.getIndex(type)][i+1] = score;
				if(DEBUG){
					System.out.println("type: " + type + "\tscore:\t" + score);
				}
				counter++;
			}
			//populated missing peaks
			double score = peakComp.compare(dummyPeak, null);
			table[this.ionTypeMap.getIndex(type)][0] = score;
			counter++;
		}
		System.out.println("Populated scoring table with total entries: " + counter);
	}
	
	private void getErrorTable(){
		int[] ranks = this.peakComp.getRankInterval();
		double[] error = this.peakComp.getMassErrorInterval();
		this.errors = new double[ranks.length][error.length];
		Peptide[] dummyPeps = {new Peptide("AAAAAAAAAA.1")};
		int counter = 0;
		for(int i = 0; i < ranks.length-1; i++){
			Peptide p = dummyPeps[0];
			MixturePeak dummyPeak = new MixturePeak(0, 0, "b", (short)1, (short)1, 0);
			dummyPeak.setPep(p);
			Peak dummyPeak2 = new Peak(0,0);
			dummyPeak2.setRank(ranks[i]);
			for(int j = 0; j < error.length-1; j++){
				//System.out.println("rank is: " + ranks[i]);
				dummyPeak2.setMoz((-1*error[j]));
				double score = this.peakComp.getErrorScore(dummyPeak, dummyPeak2);
				this.errors[i+1][j] = score;
				if(DEBUG){
					//System.out.println("rank: " + ranks[i] + "\terror" + error[j] + "\tscore:\t" + score);
					System.out.println("rank: " + i + "\terror" + error[j] + "\tscore:\t" + score);
				}
				counter++;
			}
		}
		System.out.println("Populated error scoring table with total entries: " + counter);
	}
	public double[][] getTable() {
		return table;
	}
	
	public double[][] getErrorsTable() {
		return errors;
	}
	
	public static void testAdaptor(){
		LinkedPeptidePeakScoreLearner pComp = LinkedPeptidePeakScoreLearner.loadComparator("../mixture_linked/yeast_simmix_alpha_generic_12_25.o");
		LinkedScorerAdapter adaptor = new LinkedScorerAdapter(LinkedTheoSpectrumFactory.linkedTypeMap, 
				LinkedTheoSpectrumFactory.linkedIonMap, pComp);
		adaptor.getScoringTable();
	}
	
	public static void main(String[] args){
		testAdaptor();
	}
}
