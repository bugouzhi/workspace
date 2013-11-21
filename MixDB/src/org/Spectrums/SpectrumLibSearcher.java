package org.Spectrums;
/**
 * contain various method to search a spectrum against a list of spectrum
 * @author jian wang
 *
 */
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Set;
import java.util.TreeMap;

import sequences.FastaSequence;

import mixdb.ArrayTheoreticalSpectrum;
import mixdb.MixTheoSpectrumFactory;
public class SpectrumLibSearcher {
	private SpectrumComparator filter;
	private SpectrumComparator comparator;
	private SpectrumComparator singleScorer; //only use when handling mixture spectrum where the scorer (i.e comparator) is used to score mixture
	public String spectrumFile;
	public String outputFile;
	public BufferedWriter bw;
	public static BufferedWriter defaultout=initOut();
	public SpectrumComparator getSingleScorer() {
		return singleScorer;
	}

	public void setSingleScorer(SpectrumComparator singleScorer) {
		this.singleScorer = singleScorer;
	}

	private List<Spectrum>specList;
	private List<SpectrumScorePair> spectrumScorePairs;
	
	/**
	 * Create a spectrumlib searcher with the corresponding scoring function comparator
	 * By convention, the searcher can takes two comparators, this implements many scenario
	 * where users might want to filter the whole library with a fast filter and the score
	 * the rest with a more accurate but more computationally intensive scoring function 
	 * @param specList
	 * @param comparator
	 * @param filter
	 */
	public SpectrumLibSearcher(List<Spectrum> specList, SpectrumComparator filter, SpectrumComparator comparator){
		createScorePair(specList);
		this.specList = specList;
		this.comparator = comparator;
		this.filter = filter;
		this.singleScorer  = comparator;
		if(bw == null){
			this.bw = defaultout;
		}
	}
	
	public SpectrumLibSearcher(List<Spectrum> specList, SpectrumComparator comparator){
		this(specList, comparator, comparator);
	}
	
	public static BufferedWriter initOut(){
		defaultout = new BufferedWriter(new OutputStreamWriter(System.out));
		return defaultout;
	}

	/**
	 * Given an annotated spectrum, compute number of
	 * spectrum with score better than the one with same annotation
	 * @param query
	 * @return
	 */
	
	public int rank(Spectrum query){
		sortSpecListByScore(query);
		for(int i = this.spectrumScorePairs.size()-1; i >= 0; i--){
			Spectrum curr = this.spectrumScorePairs.get(i).s;
			if(curr.getPeptide().equals(query.getPeptide())){
				return this.spectrumScorePairs.size() - i; 
			}
		}
		return -1; //only if we cannot find the correct match in the specList
	}
	
	/**
	 * similar to rank, but for mixture of a pair of peptide
	 * @param mixturequery
	 * @return
	 */
	public int[] ranks(Spectrum mixturequery){
		sortSpecListByScore(mixturequery);
		String[] peps = mixturequery.getPeptide().split(" & ");
		int rank1=-1, rank2=-1;
		for(int i = this.spectrumScorePairs.size()-1; i >= 0; i--){
			Spectrum curr = this.spectrumScorePairs.get(i).s;
			if(curr.getPeptide().equals(peps[0])){
				rank1 = this.spectrumScorePairs.size() - i; 
			}
			if(curr.getPeptide().equals(peps[1])){
				rank2 = this.spectrumScorePairs.size() - i; 
			}
			if(rank1 > 0 && rank2 > 0){
				return new int[]{rank1, rank2};
			}

		}

		return new int[]{rank1, rank2};
	}
	
	/**
	 * similar to ranks but for linked peptide
	 * @param mixturequery
	 * @return
	 */
	public int[] linkedRanks(Spectrum mixturequery){
		sortSpecListByScore(mixturequery);
		String[] peps = mixturequery.getPeptide().split("--");
		int rank1=-1, rank2=-1;
		//System.out.println("targets are: " + peps[0] + " and " + peps[1]);
		boolean found1 = false, found2 = false;
		for(int i = this.spectrumScorePairs.size()-1; i >= 0; i--){
			Spectrum curr = this.spectrumScorePairs.get(i).s;
			String pep = curr.getPeptide().split("\\.")[0];
			//System.out.println(pep + " score: " + this.spectrumScorePairs.get(i).score + " higher than target");
			if(pep.equals(peps[0]) && !found1){
				rank1 = this.spectrumScorePairs.size() - i;  
				found1 = true;
			}
			if(pep.equals(peps[1]) && !found2){
				rank2 = this.spectrumScorePairs.size() - i;  
				found2 = true;
			}
			if(rank1 > 0 && rank2 > 0){
				return new int[]{rank1, rank2};
			}


		}
		return new int[]{rank1, rank2};
	}
	
	/**
	 * search the library, find the top matching spectrum
	 * @param query
	 * @return
	 */
	public Spectrum topSpectrum(Spectrum query){
		sortSpecListByScore(query);
		if(this.spectrumScorePairs.size() <= 0){
			return null;
		}
		SpectrumScorePair best = this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1);
		//this.printTopCandidateInfo(query, best);
		System.out.println("Query " + query.peptide +  " with top answer is: " 
				+ this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1).s.peptide 
				+ " score: " + this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1).score);
		return this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1).s;
	}
	
	public Spectrum topSpectrumIter(Spectrum query){
		sortSpecListByScore(query);
		if(this.spectrumScorePairs.size() <= 0){
			return null;
		}
		SpectrumScorePair best = this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1);
		this.printTopCandidateInfo(query, best);
		TheoreticalSpectrum t = (TheoreticalSpectrum)best.s;
		removeAnnotatedPeak(query, t);
		sortSpecListByScore(query);
		best = this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1);
		this.printTopCandidateInfo(query, best);
		return this.spectrumScorePairs.get(this.spectrumScorePairs.size()-1).s;
	}
	
	
	private void removeAnnotatedPeak(Spectrum s, TheoreticalSpectrum t){
		List<Peak> toBeRemoved = getAnnotatedPeak(s, t);
		Collections.sort(toBeRemoved, PeakIntensityComparator.comparator);
		System.out.println("removing peaks: " + toBeRemoved.size() + "\t" + s.getPeak().removeAll(toBeRemoved));
		//s.getPeak().removeAll(toBeRemoved);
	}
	
	public List<Peak> getAnnotatedPeak(Spectrum s, TheoreticalSpectrum t){
		SimpleMatchingGraph g = t.getMatchGraph(s, 0.3);
		List<Peak> toBeRemoved = new ArrayList<Peak>();
		for(Iterator<Peak> it = g.vertexSet(g.Observed).iterator(); it.hasNext();){
			Peak current = it.next();
			if(g.getNeighbors(current).size() > 0){
				toBeRemoved.add(current);
			}
		}
		Collections.sort(toBeRemoved, PeakIntensityComparator.comparator);
		return toBeRemoved;
	}
	
	public Spectrum[] topSpectra(Spectrum query, int topN){
		Spectrum[] topSpectra = new Spectrum[topN];
		sortSpecListByScore(query);
		//this.printTopCandidateInfo(query, best);
		for(int i = 1; i <= topN && i <= this.spectrumScorePairs.size(); i++){
			SpectrumScorePair best = this.spectrumScorePairs.get(this.spectrumScorePairs.size()-i);
			topSpectra[i-1] = best.s;
		}
		for(int i = 1; i <= topN && i <= this.spectrumScorePairs.size(); i++){
			SpectrumScorePair best = this.spectrumScorePairs.get(this.spectrumScorePairs.size()-i);
			System.out.println("Query " + query.spectrumName + "\t" + query.peptide +  " with top answer is: " 
					+ ((TheoreticalSpectrum)this.spectrumScorePairs.get(this.spectrumScorePairs.size()-i).s).p + "\t"
					+ ((TheoreticalSpectrum)this.spectrumScorePairs.get(this.spectrumScorePairs.size()-i).s).p.getCharge()
					+ " score: " + this.spectrumScorePairs.get(this.spectrumScorePairs.size()-i).score);
		}
		
		return topSpectra;
	}
	
	
	public Spectrum[] topLinkedSpectra(Spectrum query, int topN){
		Spectrum[] topSpectra = new Spectrum[topN];
		sortSpecListByScore(query);
		int i = this.spectrumScorePairs.size()-1;
		int count = 0;
		if(this.spectrumScorePairs.size() == 0){
			return null;
		}
		double prevScore = this.spectrumScorePairs.get(i).score;
		SpectrumScorePair prev = this.spectrumScorePairs.get(i);
		while(i >= 0 && count < topN && prev.score > -100){
			SpectrumScorePair current = this.spectrumScorePairs.get(i);
			printTopLinkedCandidateInfo(query, prev);
			if(Math.abs(current.score - prevScore) > 0.000001){  //consider same score within resonable precision
				count++;
			}
			prevScore = current.score;
			prev = current;
			i--;
		}
		
		return topSpectra;
	}
	/**
	 * search the library by first use a efficient filter to filter
	 * the library, then using a more computational intensitve scorer
	 * to score the filtered library
	 * @param query
	 * @return
	 */

	public Spectrum bestSpectrum(Spectrum query){
		int topN = 1000;
		double lowScore = -100000;
		//sortSpecListByScore(query);
		int endIndex = this.spectrumScorePairs.size() - topN;
		endIndex = endIndex < 0 ? 0 : endIndex;
		List<SpectrumScorePair> filtered = new ArrayList<SpectrumScorePair>();
		System.out.println("library has size: " + this.specList.size());
		for(int i = this.spectrumScorePairs.size() - 1; i >= 0; i--){
			SpectrumScorePair curr = this.spectrumScorePairs.get(i);
			curr.s = new TheoreticalSpectrum(((TheoreticalSpectrum)curr.s).p);
			curr.score = this.comparator.compare(curr.s, query);
			filtered.add(curr);
		}
		Collections.sort(filtered);
		if(filtered.size() > 0){
			SpectrumScorePair best = filtered.get(filtered.size()-1);
			this.printTopCandidateInfo(query, best);
			return filtered.get(filtered.size()-1).s;
		}else{
			return null;
		}
	}
	
	public Spectrum bestSpectra(Spectrum query, int bestN){
		int topN = 1000;
		double lowScore = -100000;
		//sortSpecListByScore(query);
		int endIndex = this.spectrumScorePairs.size() - topN;
		endIndex = endIndex < 0 ? 0 : endIndex;
		List<SpectrumScorePair> filtered = new ArrayList<SpectrumScorePair>();
		System.out.println("library has size: " + this.specList.size());
		for(int i = this.spectrumScorePairs.size() - 1; i >= 0; i--){
			SpectrumScorePair curr = this.spectrumScorePairs.get(i);
			//curr.s = new TheoreticalSpectrum(((TheoreticalSpectrum)curr.s).p);        //why duplicated???
			curr.score = this.comparator.compare(curr.s, query);
			//System.out.println("peptide is: " + curr.s.peptide +  "\tscore:\t" + curr.score);
			filtered.add(curr);
		}
		Collections.sort(filtered);
		if(filtered.size() > 0){
			int i = filtered.size() - 1;
			while(i > 0 && i > filtered.size()- bestN){
				SpectrumScorePair best = filtered.get(i);
				this.printTopCandidateInfo(query, best);
				i--;
			}
			return filtered.get(filtered.size()-1).s;
		}else{
			return null;
		}
	}
	
	/**
	 * search the library return the best pair
	 * @param mixturequery
	 * @return
	 */
	public Spectrum[] bestPair(Spectrum mixturequery){
		sortSpecListByScore(mixturequery);
		Spectrum best = null;
		double bestScore = -1000.0, currScore = 0.0;
		int maxIndex = this.spectrumScorePairs.size()-1;
		int firstIndex = maxIndex - 5, secondIndex = maxIndex - 500;
		int besti = 0, bestj = 0;
		for(int i = maxIndex; i > firstIndex; i--){
			Spectrum s1 = this.spectrumScorePairs.get(i).s;
			for(int j = i-1; j > secondIndex; j--){
				Spectrum s2 = this.spectrumScorePairs.get(j).s;
				TheoreticalSpectrum mix = new TheoreticalSpectrum(s1.getPeptide(),  s2.getPeptide());
				currScore = this.comparator.compare(mix, mixturequery);
				best = currScore > bestScore ? mix : best;
				besti = currScore > bestScore ? i : besti;
				bestj = currScore > bestScore ? j : bestj;
				bestScore = currScore > bestScore ? currScore : bestScore;
			}
		}
		TheoreticalSpectrum single1 = new TheoreticalSpectrum(this.spectrumScorePairs.get(besti).s.getPeptide());
		TheoreticalSpectrum single2 = new TheoreticalSpectrum(this.spectrumScorePairs.get(bestj).s.getPeptide());
		double score1 = this.comparator.compare(single1, mixturequery);
		double score2 = this.comparator.compare(single2, mixturequery);
		System.out.print("Spectrum: " + mixturequery.peptide + " has best match: " + best.peptide + " with score: " 
				+ bestScore + "\t" + " : " + score1 + "\t" + score2 + "\t" + score1/single1.peptide.length() + "\t" + score2/single2.peptide.length());
		System.out.println(checkPeptidepair(best.peptide, mixturequery.peptide));
		System.out.println();
		return null;
	}
	/**
	 * search the library, return N best pairs
	 * @param mixturequery
	 * @param topN
	 * @return
	 */
	public List<Spectrum> bestCandidates(Spectrum mixturequery, int topN){
		sortSpecListByScore(mixturequery);
		//System.out.println("candidate size: " + this.spectrumScorePairs.size());
		Spectrum best = null;
		TreeMap<Double, List<Spectrum>> bestList = new TreeMap();
		double bestScore = -1000000.0, currScore = 0.0;
		int maxIndex = this.spectrumScorePairs.size()-1;
		int firstIndex = maxIndex - 20, secondIndex = maxIndex - 500;
		firstIndex = firstIndex < 0 ? 0 : firstIndex;	
		secondIndex = secondIndex < 0 ? 0 : secondIndex;
		int besti = 0, bestj = 0;
		int count=0;
		for(int i = maxIndex; i > firstIndex; i--){
			Spectrum s1 = this.spectrumScorePairs.get(i).s;
			for(int j = i-1; j >= secondIndex; j--){
				Spectrum s2 = this.spectrumScorePairs.get(j).s;
				if(this.spectrumScorePairs.get(j).score < 0 && j!= i-1){ //we make sure it gets at least one pair of canddiate
					continue;
				}
				Spectrum mix;
//				if(j == secondIndex){
					//s2 = new TheoreticalSpectrum("Z.2");
					//mix = new TheoreticalSpectrum(s1.getPeptide(),  s2.getPeptide()); //create dummy spectrum to make sure do not 
//				}else{                                                      //select a "bad" spectrum for single-peptide hits 
					mix = new TheoreticalSpectrum(s1.getPeptide(),  s2.getPeptide());
//				}
				count++;
				currScore = this.comparator.compare(mix, mixturequery);
				//System.out.println("comparing: " + mix.peptide +  "\t" + currScore);
				//currScore = Math.random();
				if(currScore > bestScore){
					List<Spectrum> cand = new ArrayList();
					cand.add(s1); 
					cand.add(s2);
					insertBestPair(cand, currScore, bestList, topN);
					if(bestList.size() >= topN){
						bestScore = bestList.firstKey().doubleValue();
					}
				}
			}
		}
		//System.out.println("best has size: " + bestList.size());
		//System.out.println("best now: " + bestList.lastKey());
		printTopCandidatesInfo(mixturequery, bestList);
		//System.out.println("total number of pairs considered: " + count);
		return null;
	}

	
	public List<Spectrum> bestArrayCandidates(Spectrum mixturequery, int topN){
		sortSpecListByScore(mixturequery);
		//System.out.println("candidate size: " + this.spectrumScorePairs.size());
		Spectrum best = null;
		TreeMap<Double, List<Spectrum>> bestList = new TreeMap();
		double bestScore = -1000000.0, currScore = 0.0;
		int maxIndex = this.spectrumScorePairs.size()-1;
		int firstIndex = maxIndex - 20, secondIndex = maxIndex - 500;
		firstIndex = firstIndex < 0 ? 0 : firstIndex;	
		secondIndex = secondIndex < 0 ? 0 : secondIndex;
		int besti = 0, bestj = 0;
		int count=0;
		for(int i = maxIndex; i > firstIndex; i--){
			Spectrum s1 = this.spectrumScorePairs.get(i).s;
			for(int j = i-1; j >= secondIndex; j--){
				Spectrum s2 = this.spectrumScorePairs.get(j).s;
				if(this.spectrumScorePairs.get(j).score < -100 && j!= i-1){ //we make sure it gets at least one pair of canddiate
					continue;
				}
				Spectrum mix;
//				if(j == secondIndex){
					//s2 = new TheoreticalSpectrum("Z.2");
					//mix = new TheoreticalSpectrum(s1.getPeptide(),  s2.getPeptide()); //create dummy spectrum to make sure do not 
//				}else{                                                      //select a "bad" spectrum for single-peptide hits 
					//mix = new TheoreticalSpectrum(s1.getPeptide(),  s2.getPeptide());
					mix = MixTheoSpectrumFactory.getMixTheoSpectrum((ArrayTheoreticalSpectrum)s1, (ArrayTheoreticalSpectrum)s2);
//				}
				count++;
				currScore = this.comparator.compare(mix, mixturequery);
				//System.out.println("comparing: " + mix.peptide +  "\t" + currScore);
				//currScore = Math.random();
				if(currScore > bestScore){
					List<Spectrum> cand = new ArrayList();
					cand.add(s1); 
					cand.add(s2);
					insertBestPair(cand, currScore, bestList, topN);
					if(bestList.size() >= topN){
						bestScore = bestList.firstKey().doubleValue();
					}
				}
			}
		}
		System.out.println("best has size: " + bestList.size());
		printTopCandidatesInfo(mixturequery, bestList);
		//System.out.println("total number of pairs considered: " + count);
		return null;
	}
	
	public List<Spectrum> bestArrayCandidates(Spectrum mixturequery, int topN, int Mix){
		sortSpecListByScore(mixturequery);
		//System.out.println("candidate size: " + this.spectrumScorePairs.size());
		Spectrum best = null;
		TreeMap<Double, List<Spectrum>> bestListPrev = new TreeMap();
		TreeMap<Double, List<Spectrum>> bestList = new TreeMap();
		double bestScore = -1000000.0, currScore = 0.0;
		int maxIndex = this.spectrumScorePairs.size()-1;
		int topK = 5, topK2 = 2000;
		int firstIndex = maxIndex - 20, secondIndex = maxIndex - 500;
		int count = 0;
		firstIndex = firstIndex < 0 ? 0 : firstIndex;	
		secondIndex = secondIndex < 0 ? 0 : secondIndex;
		int besti = 0, bestj = 0;
		for(int i = 0; i < topK; i++){
			List<Spectrum> cand = new ArrayList();
			cand.add(this.spectrumScorePairs.get(maxIndex-i).s);
			SpectrumScorePair curr = this.spectrumScorePairs.get(maxIndex-i);
			insertBestPair(cand, curr.score, bestListPrev, topK);
		}

		for(int m = 2; m <= Mix; m++){
			bestList = new TreeMap();
			bestScore = -10000;
			for(Iterator it = bestListPrev.keySet().iterator(); it.hasNext();){
				Spectrum s = bestListPrev.get(it.next()).get(0);
				for(int j = maxIndex; j >= maxIndex - topK2; j--){
					Spectrum s2 = this.spectrumScorePairs.get(j).s;
					Spectrum mix = MixTheoSpectrumFactory.getMixTheoSpectrum((ArrayTheoreticalSpectrum)s, (ArrayTheoreticalSpectrum)s2, m);
					count++;
					currScore = this.comparator.compare(mix, mixturequery);
					if(currScore > bestScore){
						List<Spectrum> cand = new ArrayList();
						cand.add(mix);
						insertBestPair(cand, currScore, bestList, topK);
						if(bestList.size() >= topK){
							bestScore = bestList.firstKey().doubleValue();
						}
					}
				}
			    //System.out.println("total combination so-far at this round: " + count);
			}
			//System.out.println("best-" + m + " has size: " + bestList.size());
			bestListPrev = bestList;
		}
		System.out.println("total combiniations considered: " + count);
		System.out.println("best has size: " + bestList.size());
		printMultiplexTopCandidatesInfo(mixturequery, bestList);
		//System.out.println("total number of pairs considered: " + count);
		return null;
	}
	

	
	private void printTopCandidateInfo(Spectrum query, SpectrumScorePair match){
		//TheoreticalSpectrum th = (TheoreticalSpectrum)match.s;
		TheoreticalSpectrum th = new TheoreticalSpectrum(match.s.getPeptide());
		double[] stat = th.analyzeAnnotation(query, match.s.peptide, 0.3, false);
		if(match.s instanceof ArrayTheoreticalSpectrum){
			ArrayTheoreticalSpectrum arry = (ArrayTheoreticalSpectrum)match.s;
			FastaSequence seq = arry.getPeplite().getFastaseq();
			String annot = seq.getAnnotation(arry.getPeplite().getBeginInd()).split("\\s+")[0];
			System.out.println(this.spectrumFile + "\t" + query.scanNumber + "\t" +  th.p + "\t" + annot   
					+"\t" + query.parentMass + "\t"  + query.charge + "\t" + match.s.parentMass + "\t" + th.charge + "\t" + match.score + "\t" 
					+ (match.score/match.s.peptide.length()) + "\t" + stat[0] + "\t" + stat[1] + "\t" + stat[2] + "\t" 
					+ stat[3] + "\t" + stat[4] + "\t" + stat[5]);
		}else{
			System.out.println(this.spectrumFile +"\t" + query.getPeptide() +  "\t" +  th.p 
			+"\t" + query.parentMass + "\t"  + query.charge + "\t" + match.s.parentMass + "\t" + th.charge + "\t" + match.score + "\t" 
			+ (match.score/match.s.peptide.length()) + "\t" + stat[0] + "\t" + stat[1] + "\t" + stat[2] + "\t" 
			+ stat[3] + "\t" + stat[4] + "\t" + stat[5]);
		}
	}
	
	private void printTopCandidatesInfo(Spectrum query, TreeMap<Double, List<Spectrum>> bestList){
		for(Iterator<Double> it = bestList.keySet().iterator(); it.hasNext();){
			Double key = it.next();
			List<Spectrum> best = bestList.get(key);
			double score1 = this.singleScorer.compare(best.get(0), query);
			double score2 = this.singleScorer.compare(best.get(1), query);
			//TheoreticalSpectrum th1 = new TheoreticalSpectrum(best.get(0).peptide, new String[]{"b"}, new String[]{"y"});
			//TheoreticalSpectrum th2 = new TheoreticalSpectrum(best.get(1).peptide, new String[]{"b"}, new String[]{"y"});
			//double score1 = this.singleScorer.compare(th1, query);
			//double score2 = this.singleScorer.compare(th2, query);
			TheoreticalSpectrum th = new TheoreticalSpectrum(best.get(0).peptide, best.get(1).peptide);
			String p1 = best.get(0).peptide.split("\\.")[0];
			String p2 = best.get(1).peptide.split("\\.")[0];
			if(p1.equals(p2)){ //we do not allow same peptides 
				//continue;
			}
			double[] stat = th.analyzeMixtureAnnotation(query, best.get(0).peptide, best.get(1).peptide, 0.3);
			String bestpeptide = best.get(0).peptide + "\t" + best.get(1).peptide;
			try{
			if(best.get(0) instanceof ArrayTheoreticalSpectrum){
				ArrayTheoreticalSpectrum arry1 = (ArrayTheoreticalSpectrum)best.get(0);
				ArrayTheoreticalSpectrum arry2 = (ArrayTheoreticalSpectrum)best.get(1);
				FastaSequence seq = arry1.getPeplite().getFastaseq();
				String annot1 = seq.getAnnotation(arry1.getPeplite().getBeginInd()).split("\\s+")[0];
				String annot2 = seq.getAnnotation(arry2.getPeplite().getBeginInd()).split("\\s+")[0];
				th = new TheoreticalSpectrum(arry1.getPeplite().toString() +"." + arry1.charge, 
						arry2.getPeplite().toString()+"."+arry2.charge);
				stat = th.analyzeMixtureAnnotation(query, arry1.getPeplite().toString(), arry2.getPeplite().toString(), 0.3);
				bestpeptide = arry1.getPeplite().toString() + "\t" + arry2.getPeplite().toString();
				String pep ="";
				if(!query.peptide.contains("Dummy")){
					pep = "\t"+query.peptide;
				}
				bw.write(this.spectrumFile + "\t" + query.scanNumber + "\t" + query.parentMass + "\t" + best.get(0).parentMass + "\t"  + best.get(1).parentMass + "\t" + query.charge
						+ "\t" +  + arry1.charge + "\t"  + arry2.charge + "\t" + bestpeptide  + "\t"  + annot1 + "\t" + annot2 
						+ " \t" + key.doubleValue() + "\t"  + score1 + "\t" + score2 + "\t" + score1/best.get(0).peptide.length() + "\t" + score2/best.get(1).peptide.length()
						+ "\t" + stat[0] + "\t" + stat[1] + "\t"
						+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
						+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
						+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12]+ "\n");
			}else if(query.peptide.contains(" & ")){
				bw.write("Spectrum:\t" + query.scanNumber + "\t" + query.getPeptide() + " best: " +  query.parentMass + "\t" + best.get(0).parentMass + "\t"  + best.get(1).parentMass 
					+ "\t" +  bestpeptide 
					+ " \t" + key.doubleValue() + "\t"  + score1 + "\t" + score2 + "\t" + score1/best.get(0).peptide.length() + "\t" + score2/best.get(1).peptide.length()
					+ "\t" + stat[0] + "\t" + stat[1] + "\t"
					+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
					+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
					+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12] + "\n");
			}else{
				System.out.print("Spectrum:\t" + query.scanNumber + "\t" + " best:\t" +  query.parentMass + "\t" + best.get(0).parentMass + "\t"  + best.get(1).parentMass 
						+ "\t" +  bestpeptide  + " \t" 
						+ key.doubleValue() + "\t"  + score1 + "\t" + score2 + "\t" + score1/best.get(0).peptide.length() + "\t" + score2/best.get(1).peptide.length()
						+ "\t" + stat[0] + "\t" + stat[1] + "\t"
						+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
						+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
						+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12] + "\n");
			}
			bw.flush();
			}catch(IOException ioe){
				System.err.println(ioe.getMessage());
				ioe.printStackTrace();
			}
		}
		//System.out.println();
	}
	
	private void printMultiplexTopCandidatesInfo(Spectrum query, TreeMap<Double, List<Spectrum>> bestList){
		for(Iterator<Double> it = bestList.keySet().iterator(); it.hasNext();){
			Double key = it.next();
			List<Spectrum> best = bestList.get(key);
			try{
				bw.write(this.spectrumFile + "\t" + query.scanNumber + "\t" + query.parentMass + "\t" + query.charge +"\t" + key + "\t" + best.get(0).peptide);
				bw.write("\n");
				bw.flush();
			}catch(IOException ioe){
				System.err.println(ioe.getMessage());
				ioe.printStackTrace();
			}
		}
		//System.out.println();
	}
	
	
	private void printTopLinkedCandidateInfo(Spectrum query, SpectrumScorePair match){
			TheoreticalSpectrum th = (TheoreticalSpectrum)match.s;
			if(th instanceof LazyEvaluateLinkedSpectrum){
				((LazyEvaluateLinkedSpectrum) th).createSpectrum();
			}
			String[] peptides = th.peptide.split(" & ");
			LinkedPeptide lp = (LinkedPeptide)th.p;
			TheoreticalSpectrum t1 = new TheoreticalSpectrum(lp.peptides[0], lp.peptides[0].getCharge());
			TheoreticalSpectrum t2 = new TheoreticalSpectrum(lp.peptides[1], lp.peptides[1].getCharge());
			double score1 = this.singleScorer.compare(t1, query);
			double score2 = this.singleScorer.compare(t2, query);
			if(t1.peptide.equals(t2.peptide)){ //we do not allow same peptides 
				return;
			}
			double[] stat = th.analyzeMixtureAnnotation(query, t1.peptide, t2.peptide);
			String bestpeptide = th.peptide;
			String peptide1 = peptides[0].replaceAll("[\\.\\+0-9]", "");
			String peptide2 = peptides[1].replaceAll("[\\.\\+0-9]", "");
			System.out.print("Spectrum: " + query.spectrumName + "\tbest:\t" + query.getPeptide() + "\t" + query.parentMass +"\t" + query.charge +
					"\t"+  bestpeptide + "\t" + match.s.parentMass + "\t" + match.score + "\t"  + score1 + "\t" + score2 + "\t" + score1/peptide1.length() + "\t" + score2/peptide2.length()
					+ "\t" + stat[0] + "\t" + stat[1] + "\t"
					+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
					//+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
					+ stat[13] + "\t" + stat[14] + "\t" + stat[15] + "\t" + stat[16] + "\t"
					+ stat[9] + "\t" + stat[10] + "\t" + stat[11] + "\t" + stat[12]);
			//System.out.println("\t" + checkPeptidepair(bestpeptide, query.peptide));	
			System.out.println();
	}
	
	private void insertBestPair(List<Spectrum> candidate, double score, TreeMap<Double, List<Spectrum>> bestPairs, int numKept){
		if(bestPairs.containsKey(score)){
			//score = score + 0.00000001;  //avoid duplicate scores   @TODO maybe should add check to make sure the canddiate pairs are not the same
		}
		
		if(bestPairs.size() >= numKept+1){
			bestPairs.pollFirstEntry();
		}
		bestPairs.put(new Double(score), candidate);
		//System.out.println("best now: " + bestPairs.lastKey());
	}
	
	public Spectrum[] bestLinkedPair(Spectrum mixturequery){
		sortSpecListByScore(mixturequery);
		Spectrum best = null;
		double bestScore = -1000.0, currScore = 0.0;
		int maxIndex = this.spectrumScorePairs.size()-1;
		int firstIndex = maxIndex - 20, secondIndex = maxIndex - 600;
		int besti = 0, bestj = 0;
		if(this.specList.size() < 2){
			return null;
		}
		for(int i = maxIndex; i > (firstIndex > 0 ? firstIndex : -1); i--){
			TheoreticalSpectrum s1 = (TheoreticalSpectrum)this.spectrumScorePairs.get(i).s;
			for(int j = i-1; j > (secondIndex > 0 ? secondIndex : -1); j--){
				TheoreticalSpectrum s2 = (TheoreticalSpectrum)this.spectrumScorePairs.get(j).s;
				if(Math.abs(s1.parentMass - s2.parentMass) < 0.05){							
					TheoreticalSpectrum mix = new TheoreticalSpectrum(s1.p,  s2.p, (short)mixturequery.charge, true);
					currScore = this.comparator.compare(mix, mixturequery);
					//System.out.println(mix.getPeptide() + " score: " + currScore);
					best = currScore > bestScore ? mix : best;
					besti = currScore > bestScore ? i : besti;
					bestj = currScore > bestScore ? j : bestj;
					bestScore = currScore > bestScore ? currScore : bestScore;
				}
			}
		}
		TheoreticalSpectrum single1 = new TheoreticalSpectrum(this.spectrumScorePairs.get(besti).s.getPeptide());
		TheoreticalSpectrum single2 = new TheoreticalSpectrum(this.spectrumScorePairs.get(bestj).s.getPeptide());
		double score1 = this.comparator.compare(single1, mixturequery);
		double score2 = this.comparator.compare(single2, mixturequery);
		System.out.print("Spectrum: " + mixturequery.peptide + " has best match: " + best.peptide + " with score:\t" + bestScore + "\t" + " : " + score1 + "\t" + score2);
		//System.out.println(checkPeptidepair(best.peptide, mixturequery.peptide));
		System.out.println();
		return null;
	}
	
	public boolean checkPeptidepair(String peptide1, String peptide2){
		String[] peps1 = peptide1.split(" & ");
		String[] peps2 = peptide2.split(" & ");
		return (peps1[0].equals(peps2[0]) && peps1[1].equals(peps2[1]))
			|| (peps1[0].equals(peps2[1]) && peps1[1].equals(peps2[0]));
	}
	
	private void sortSpecListByScore(Spectrum query){
		//System.out.println("Scan: " + query.scanNumber +  "\tlibrary has size: " + this.specList.size());
		for(Iterator<SpectrumScorePair> it = spectrumScorePairs.iterator(); it.hasNext();){
			SpectrumScorePair curr = it.next();
			curr.score = this.filter.compare(curr.s, query);
			//System.out.println(((ArrayTheoreticalSpectrum)curr.s).getPeplite() +  "\t" + curr.s.parentMass + "\t score is: " + curr.score);
		}
		Collections.sort(this.spectrumScorePairs);
	}
	
	//note: we can avoid excessive temp obj creation by using a static list
	//to old and recycle spectrumScore pair
	private void createScorePair(List<Spectrum> specList){
		this.spectrumScorePairs = new ArrayList(specList.size());
		for(Iterator<Spectrum> it = specList.iterator(); it.hasNext();){
			Spectrum currentSpect = it.next();
			this.spectrumScorePairs.add(new SpectrumScorePair(currentSpect));
		}
	}
	
	//a temporary container holding a score for each spectrum in the list
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
}
