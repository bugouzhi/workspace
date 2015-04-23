package mixdb;

import java.util.Iterator;
import java.util.Map;

import org.Spectrums.LabelledPeak;
import org.Spectrums.Peak;
import org.Spectrums.Peptide;
import org.Spectrums.RankBaseScoreLearner;

/**
 * Convert old scoring function to new one (see TheoreticalSpectrumFactrory)
 * This adaptor provide methods that covert scoring models in the old framework
 * to the new one.
 * @author Jian Wang
 *
 */
public class RankBaseScoreAdapter {
	protected Map<String, IonType> typeMap;
	protected IonTypeMapper ionTypeMap;
	protected RankBaseScoreLearner peakComp;
	protected double[][] table; //the new scoring table
	protected double[][] errors; //mass error scores
	protected boolean DEBUG=false;
	
	public RankBaseScoreAdapter(Map<String, IonType> typeMap, IonTypeMapper ionMap,RankBaseScoreLearner pComp){
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
		for(Iterator<IonType> it = typeMap.values().iterator(); it.hasNext();){
			SimpleIonType type = (SimpleIonType)it.next();
			SimplePeptideType pType = (SimplePeptideType)type.getPType(); 
			Peptide p = dummyPeps[pType.getLength()];
			p.setCharge((short)pType.getPepCharge());
			LabelledPeak dummyPeak = new LabelledPeak(0, 0, type.getType(), (short)1, (short)type.getCharge());
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
			LabelledPeak dummyPeak = new LabelledPeak(0, 0, "b", (short)1, (short)2);
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
		RankBaseScoreLearner pComp = RankBaseScoreLearner.loadComparator("../mixture_linked/yeast_single_model_realannotated_win10_25.o");
		RankBaseScoreAdapter adaptor = new RankBaseScoreAdapter(TheoreticalSpectrumFactory.standardTypeMap, 
				TheoreticalSpectrumFactory.standardIonMap, pComp);
		//adaptor.getScoringTable();
	}
	
	public static void main(String[] args){
		testAdaptor();
	}
	
	

}
