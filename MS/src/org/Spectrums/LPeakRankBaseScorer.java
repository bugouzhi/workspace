package org.Spectrums;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import Utils.ArrayUtils;

public class LPeakRankBaseScorer implements PeakComparator, Serializable{
	private static final long serialVersionUID = 193712381028301L;
	private RankBaseScoreLearner core;
	public LPeakRankBaseScorer(SpectrumLib lib){
		this.core = new RankBaseScoreLearner(lib);
		this.core.getIonsCount();
	}
	
	public LPeakRankBaseScorer(String libFile){
		this.core = new RankBaseScoreLearner(libFile);
		this.core.getMixtureIonCount();
	}
	
//	public LPeakScoreLearner(SpectrumLib lib, String[] ionsType, 
//			int minCharge, int maxCharge){
//		super(lib, ionsType, minCharge, maxCharge);
//		getIonsCount();
//	}
	
	private int[] getLinkedIndex(LabelledPeak lp, Peak realPeak){
		int rankIndex;
		if(realPeak == null){
			rankIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.core.getRankInterval())+1;
		}
		if(lp == null){
			return new int[]{0,0,0,0,rankIndex};
		}
		Peptide p = lp.getPep();
		
		
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.core.getLengthInterval());
		///System.out.println("peptide is: " + p.getPeptide() + " index : " + peptideLength);
		int peptideCharge = getPepCharge(lp)-1;
		int peakCharge = getPeakCharge(lp)-1;
		int ionIndex = this.core.getIonIndex(lp)+1;
		return new int[]{peptideLength, peptideCharge, peakCharge, ionIndex, rankIndex,0};
	}
	
	
	private int getPepCharge(LabelledPeak lp){
		Peptide pep = lp.getPep();
		return getPepCharge(pep);
	}
	
	private int getPepCharge(Peptide p){
//		return p.getCharge();
		if(p.getCharge() < 4){
			return 2;
		}else{
			return 3;
		}
	}
	
	private int getPeakCharge(LabelledPeak lp){
		Peptide pep = lp.getPep();
		if(TheoreticalSpectrum.isLinkedPeak(pep, lp)){
			if(pep.getCharge() == 2){
				return lp.getCharge();
			}else if(pep.getCharge() == 3){
				return lp.getCharge() - 1;
			}else if(pep.getCharge() == 4){
				return lp.getCharge() - 1;
			}else if(pep.getCharge() == 5){
				return lp.getCharge() - 2;
			}else{
				if(lp.getCharge() > 5){
					return lp.getCharge()-3;
				}else{
					return lp.getCharge() - 2;
				}
			}
		}else{
			return lp.getCharge();
		}
	}
	
	
	private int[] getLinkedNoiseIndex(Peptide p, Peak realPeak){
		int rankIndex;
		if(realPeak == null){
			rankIndex = 0;
		}else{
			rankIndex = ArrayUtils.getIntervalIndex(realPeak.getRank(), this.core.getRankInterval())+1;
		}
		int peptideLength = ArrayUtils.getIntervalIndex(p.getPeptide().length(), this.core.getLengthInterval());
		int peptideCharge = getPepCharge(p)-1;
		int peakCharge = 0;
		int ionIndex = 0;
		return new int[]{peptideLength, peptideCharge, peakCharge, ionIndex, rankIndex,0};
	}
	
	@Override
	public double compare(Peak p1, Peak p2) {
		LabelledPeak lp = (LabelledPeak)p1;
		if(p1 == null ){
			return 0;
		}else{
			int[] index =  getLinkedIndex(lp, p2);
			int[] index2 = getLinkedNoiseIndex(lp.getPep(), p2); //assuem if we were to match the peak to noise
			double score = this.core.getValue(index);
			double score2 = this.core.getValue(index2);
			int[] errorIndex = this.core.getErrorIndex(lp, p2);
			if(p2 != null){
				double errorScore = this.core.getErrorModel().get(errorIndex);
				score *= errorScore;
				if(score == 0){
					//score = 0.00001;
					return 0;
				}
				score2 = score2 / (this.core.getMassErrorInterval().length-1);
				//System.out.println("peptide " + "\t" +  lp +  "\trank: " + p2+ "\t" + p2.getRank() + " score: " + score + "\t" + score2 + "\t" + Math.log(score/score2) + "\t" + errorScore);
				if(score2 == 0){
					//System.out.println("warning noise is zero: " + lp + "\t" + p2);
					return 0;
				}

			}
//			System.out.println("score: " + score);
//			System.out.println("score2: " + score2);
			if(Double.isNaN(score)){
				return 0;
			}
			return Math.log(score/score2);
		}
	}
	
	public void writeLibToFile(String outfile){
		try{
			BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(outfile));
			ObjectOutputStream oo = new ObjectOutputStream(bo);
		    oo.writeObject(this);
		    oo.flush();
		    oo.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		}
	}
		
	public static LPeakRankBaseScorer loadComparator(String file){
		try{
			BufferedInputStream bi = new BufferedInputStream(new FileInputStream(file));
			ObjectInputStream oi = new ObjectInputStream(bi);
		    Object o = oi.readObject();
		    return (LPeakRankBaseScorer)o;
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
}
