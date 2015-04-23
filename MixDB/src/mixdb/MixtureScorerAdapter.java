package mixdb;

import java.util.Iterator;
import java.util.Map;

import org.Spectrums.LabelledPeak;
import org.Spectrums.MixturePeak;
import org.Spectrums.MixturePeakScoreLearner;
import org.Spectrums.Peak;
import org.Spectrums.Peptide;
import org.Spectrums.RankBaseScoreLearner;
import org.Spectrums.Spectrum;
/**
 * Please see TheoreticalSpectrumFactory for the basic framework of theoretical spectrum.
 * This adaptor class provide methods to transform the scoring model under the old framework to our new framework so 
 * scoring can be done under our new framework using learned scoring models.
 * @TODO
 * Since the new framework for theoretical spectrum is much cleaner.  It wil be greate if we
 * add in functionality that can learn the scoring mode only base on this new framework  
 * @author Jian
 *
 */
public class MixtureScorerAdapter{
	private Map<String, MixIonType> typeMap;
	private IonTypeMapper ionTypeMap;
	private MixturePeakScoreLearner peakComp;
	private double[][] table; //the new scoring table
	private double[][] errors;
	private boolean DEBUG=false;
	
	public MixtureScorerAdapter(Map<String, MixIonType> typeMap, IonTypeMapper ionMap, MixturePeakScoreLearner pComp){
		this.typeMap = typeMap;
		this.ionTypeMap = ionMap;
		this.peakComp = pComp;
		getScoringTable();
		getErrorTable();
	}
	
	private void getScoringTable(){
		int[] ranks = this.peakComp.getRankInterval();
		this.table = new double[typeMap.values().size()][ranks.length];
		Peptide[] dummyPeps = {new Peptide("AAAAAAAAA.1"), new Peptide("AAAAAAAAAAAAA.1")};
		int counter = 0;
		for(Iterator<MixIonType> it = typeMap.values().iterator(); it.hasNext();){
			MixIonType type = (MixIonType)it.next();
			Spectrum s = new Spectrum();
			s.charge = type.getMixtureCharge();
			MixturePeptideType pType = (MixturePeptideType)type.getPType(); 
			Peptide p = dummyPeps[pType.getLength()];
			p.setCharge((short)pType.getPepCharge());
			MixturePeak dummyPeak = new MixturePeak(0, 0, type.getType(), (short)1, (short)type.getCharge(), pType.getPeptideAbundance());
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
		Peptide[] dummyPeps = {new Peptide("AAAAAAAAA.1"), new Peptide("AAAAAAAAAAAAA.1")};
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
		MixturePeakScoreLearner pComp = MixturePeakScoreLearner.loadComparator("../mixture_linked/yeast_simmix_alpha_generic_12_25.o");
		MixtureScorerAdapter adaptor = new MixtureScorerAdapter(MixTheoSpectrumFactory.mixTypeMap, 
				MixTheoSpectrumFactory.mixIonMap, pComp);
		adaptor.getScoringTable();
	}
	
	public static void main(String[] args){
		testAdaptor();
	}
}
