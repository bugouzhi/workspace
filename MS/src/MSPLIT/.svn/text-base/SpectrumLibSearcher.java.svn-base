package MSPLIT;
/**
 * contain various method to search a spectrum against a list of spectrum
 * @author jian wang
 *
 */
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;


public class SpectrumLibSearcher {
	private SpectrumComparator filter;
	private SpectrumComparator comparator;
	private List<Spectrum>specList;
	private List<SpectrumScorePair> spectrumScorePairs;
	private SpectrumScorePairFactory factory;
	
	/**
	 * Create a spectrumlib searcher with the corresponding scoring function comparator
	 * By convention, the searcher can takes two comparators, this implements many scenario
	 * where users might want to filter the whole library with a fast filter and the score
	 * the rest with a more accurate but more computationally intensive scoring function 
	 * @param specList
	 * @param comparator
	 * @param filter
	 */

	public SpectrumLibSearcher(List<Spectrum> specList, SpectrumComparator comparator){
		this.filter = comparator;
		this.comparator = comparator;
		this.specList = specList;
		this.factory = new SpectrumScorePairFactory();
		this.createScorePair(this.specList);
	}
	
	public List<SpectrumScorePair> topCandidates(Spectrum query, int topN){
		sortSpecListByScore(query);
		int toIndex = this.spectrumScorePairs.size();
		int fromIndex = toIndex - topN;
		fromIndex = fromIndex > 0 ? fromIndex : 0;
		return this.spectrumScorePairs.subList(fromIndex, toIndex);
	}
	
	public SpectrumScorePair topCandidate(Spectrum query){
		sortSpecListByScore(query);
		SpectrumScorePair best = this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1);
		printTopCandidateInfo(query, best);
		return best;
	}
	
	public void printTopCandidateInfo(Spectrum query, SpectrumScorePair top){
		System.out.println(query.getSpectrumName() + "\t" 
				+ query.getParentMass() + "\t" 
				+ query.getCharge() + "\t"
				+ top.score + "\t"
				+  ((SimpleAnnotatedSpectrum)top.s).getPeptide() + "\t"
				+ top.s.getParentMass() + "\t" 
				+ top.s.getCharge());
	}
	
	public void printTopCandidatePairInfo(Spectrum query, SpectrumScorePair top){
		AnnotatedMixtureSpectrum best = (AnnotatedMixtureSpectrum)top.s;
		System.out.println(query.getSpectrumName() + "\t" 
				+ query.getParentMass() + "\t" 
				+ query.getCharge() + "\t"
				+ top.score + "\t" 
				+ best.getPeptides()[0] + "\t" + best.getPeptides()[1] + "\t"
				+ best.getParentMasses()[0] + "\t" + best.getParentMasses()[1] + "\t"
				+ best.getCharges()[0] +"\t" + best.getCharges()[1]);
	}
	//get the top-scoring pair, rigourously speaking, topN pairs
	//may not be available since a branch and bound strategies is used
	//to look for the best pairs, only topN-scoring pairs encounter
	//during the search is possible for return, and not sure if this
	//acutally means something
	private double getUpperBound(SpectrumScorePair p1, SpectrumScorePair p2){
		double S1 = p1.score*p1.score;
		double S2 = p2.score*p2.score;
		return (S1 + S2)/Math.pow(S1, S2);
	}
	
	public Spectrum topCandidatePair(Spectrum query){
		sortSpecListByScore(query);
		SpectrumScorePair curr1, curr2, top, best1, best2;
		Spectrum curr, best;
		Spectrum original = query; //save the original mixture spectrum
		// mix = mix.toNormVector(1, 0.5, 2000); //we make sure we normalize the mixture spectrum, since the bound makes this assumption
		double bestscore=0, nextscore = 0, upperbound = 0;
		double alpha = 1, S1=0, S2=0;
		double bestalpha = 0;
		int counts = 0; //we keep track of combo we need to create in order to track the efficacy of our bound
		top = this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1);
		List<SpectrumScorePair> candidates = new ArrayList();
		candidates.addAll(this.spectrumScorePairs);
		best1 = null; best = null; best2 = null;
		for(int i = candidates.size()-1; i > 0; i--){
			//System.out.println("counts: " + counts);
			curr1 = candidates.get(i);
			AnnotatedSpectrum s1 = (AnnotatedSpectrum)curr1.s;
			if(counts > 50000){
				break;
			}
			for(int j = i-1; j > 0; j--){
				curr2 = candidates.get(j);
				AnnotatedSpectrum s2 = (AnnotatedSpectrum) curr2.s;
				if(counts > 50000){
					break;
				}
				if(getUpperBound(curr1, curr2) > bestscore){
					alpha = CosineSpectrumComparator.alpha(query, curr1.s, curr2.s);
					//System.out.println("alpha is: " + alpha);
					Spectrum mix1 = new AnnotatedMixtureSpectrum(s1, s2, 1, alpha);
					Spectrum mix2 = new AnnotatedMixtureSpectrum(s1, s2, 1, alpha);
					double score1 = this.comparator.compare(query, mix1);
					double score2 = this.comparator.compare(query, mix2);
					curr = score1 > score2 ? mix1 : mix2;
					nextscore = score1 > score2 ? score1 : score2;
					//System.out.println("current score: " + nextscore);
					if(nextscore > bestscore){
						bestscore = nextscore;
						best1 = candidates.get(i);
						best2 = candidates.get(j);
						bestalpha = alpha;
						best = curr;
					}
					//System.out.println("i: " + i + " j: " + j);
					counts++;
				}else{
					break;
				}
			}
		}
		if(bestalpha < 1){
			bestalpha = 1/bestalpha;
		}
		SpectrumScorePair bestCand = new SpectrumScorePair(best);
		bestCand.score = bestscore;
		printTopCandidatePairInfo(query, bestCand);
		return best;
	}

	public List<SpectrumScorePair> topCandidatePair(Spectrum query, int topN){
		return null;
	}

	//obtain statistics for top match candidate to calculate significance
	//of the match
	public double[] getTopMatchStat(Spectrum query, Spectrum bestPair, 
				Spectrum best1, Spectrum best2, double bestScore, double alpha){
		double score1 = this.comparator.compare(query, best1);
		double score2 = this.comparator.compare(query, best2);
		double pscore1 = CosineSpectrumComparator.ProjectedCosineComparator.compare(query, best1);
		double pscore2 = CosineSpectrumComparator.ProjectedCosineComparator.compare(query, best2);
		double resAlpha = SpectrumUtil.residual(query, best1);
		int intensePeak = SpectrumUtil.explainedIntensity(query, 0.85);
		return null;
	}
	
	private void sortSpecListByScore(Spectrum query){
		System.out.println("library has size: " + this.specList.size());
		for(Iterator<SpectrumScorePair> it = spectrumScorePairs.iterator(); it.hasNext();){
			SpectrumScorePair curr = it.next();
			curr.score = this.comparator.compare(curr.s, query);
			//System.out.println(" score is: " + curr.score);
		}
		Collections.sort(this.spectrumScorePairs);
	}

	private void createScorePair(List<Spectrum> specList){
		this.spectrumScorePairs = new ArrayList(specList.size());
		for(Iterator<Spectrum> it = specList.iterator(); it.hasNext();){
			Spectrum currentSpect = it.next();
			this.spectrumScorePairs.add(new SpectrumScorePair(currentSpect));
		}
	}
	
	
	/**
	 * a temporary container holding a score for each spectrum in the list
	 */
	class SpectrumScorePair implements Comparable<SpectrumScorePair>{
		public Spectrum s;
		public double score;
		public SpectrumScorePair(Spectrum s){
			this.s = s;
		}
		@Override
		public int compareTo(SpectrumScorePair s) {
			double diff= this.score - s.score;
			if(diff > 0){
				return 1;
			}else if(diff == 0){
				return 0;
			}else{
				return -1;
			}
		}
	}
	
	/**
	 * A factory for SpectrumScore pair, can re-cycle these temporary objects
	 * to increase efficiency
	 * @author bugouzhi
	 *
	 */
	class SpectrumScorePairFactory{
		private List<SpectrumScorePair> pairList; 
		int currentIndex=0;
		public SpectrumScorePairFactory(){
			this.pairList = new ArrayList<SpectrumScorePair>();
		}
		public SpectrumScorePair createScorePair(Spectrum s){
			if(currentIndex > pairList.size()){
				SpectrumScorePair p = new SpectrumScorePair(s);
				pairList.add(p);
				currentIndex++;
				return p;
			}else{
				return pairList.get(currentIndex);
			}
		}
		
		public void resetFactory(){
			this.currentIndex = 0;
		}
		
	}
}
