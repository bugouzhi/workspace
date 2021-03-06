package msgf;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class MSGFDBResultGenerator extends ArrayList<MSGFDBResultGenerator.DBMatch> {
	/**
	 * 
	 */
	
	private static final int NUM_SPECS_TO_USE_SIMPLE_ETDA_FORMULA = 30000;
	private static final long serialVersionUID = 1L;
	private String header;
	public MSGFDBResultGenerator(String header)
	{
		this.header = header;
	}
	public void computeEFDR()
	{
		double cumulativePValue = 0;
		boolean useComplicatedFormula = true;
		if(this.size() >= NUM_SPECS_TO_USE_SIMPLE_ETDA_FORMULA)
			useComplicatedFormula = false;
		for(int i=0; i<this.size(); i++)
		{
			double specProb = get(i).getSpecProb();
			double pValue = get(i).getPValue();
			cumulativePValue += pValue;
			double eTD = (i+1) - cumulativePValue;	// expected target discovery
			double eDD = cumulativePValue;	// expected decoy discovery
			if(useComplicatedFormula)
			{
				for(int j=i+1; j<this.size(); j++)
					eDD += get(j).getEDD(specProb);
			}
			else
			{
				eDD += pValue*(this.size()-(i+1));
			}
			get(i).setEFDR(Math.min(eDD/eTD, 1));
		}
	}
	
	public void writeResults(PrintStream out, boolean printEFDR)
	{
		if(printEFDR)
			out.println(header+"\tEFDR");
		else
			out.println(header);
		String eFDRStr;
		for(MSGFDBResultGenerator.DBMatch m : this)
		{
			if(printEFDR)
			{
				double eFDR = m.getEFDR();
				if(eFDR < Float.MIN_NORMAL)
					eFDRStr = String.valueOf(eFDR);
				else
					eFDRStr = String.valueOf((float)eFDR);
				out.println(m.getResultStr()+"\t"+eFDRStr);
			}
			else
				out.println(m.getResultStr());
		}
	}
	
	public static class DBMatch implements Comparable<DBMatch>
	{
		private double specProb;
		private double pValue;
		private int numPeptides;
		private String resultStr;
		private double[] cumScoreDist;
		private double eFDR;
		int curIndex;
		
		public DBMatch(double specProb, int numPeptides, String resultStr, ScoreDist scoreDist) {
			this.specProb = specProb;
			this.pValue = getPValue(specProb, numPeptides);
			this.numPeptides = numPeptides;
			this.resultStr = resultStr;
			
			if(scoreDist != null && scoreDist.isProbSet())
			{
				this.cumScoreDist = new double[scoreDist.getMaxScore()-scoreDist.getMinScore()+1];
				cumScoreDist[0] = 0;
				int index = 1;
				for(int t=scoreDist.getMaxScore()-1; t>=scoreDist.getMinScore(); t--)
				{
					cumScoreDist[index] = cumScoreDist[index-1] + scoreDist.getProbability(t);
					index++;
				}
			}
			curIndex = 0;
		}
		
		public static double getPValue(double specProb, int numPeptides)
		{
			double pValue;
			double probCorr = 1.-specProb;
			if(probCorr < 1.)
				pValue = 1.- Math.pow(probCorr, numPeptides);
			else
				pValue = specProb*numPeptides;
			return pValue;
		}
		
		public void setEFDR(double eFDR)	{ this.eFDR = eFDR; }
		
		public double getEFDR() {
			return eFDR;
		}
		
		/**
		 * Gets expected decoy discovery for a given specProbThreshold
		 */
		public double getEDD(double specProbThreshold)	{
			double probEqualOrBetterTargetPep;
			if(specProbThreshold >= specProb)
				probEqualOrBetterTargetPep = specProb;
			else
				probEqualOrBetterTargetPep = getSpectralProbability(specProbThreshold);
			
			double pValue = getPValue(probEqualOrBetterTargetPep, numPeptides);
			return pValue;
		}
		
		// returns cumulative probability <= specProbThreshold
		public double getSpectralProbability(double specProbThreshold)
		{
//			int index = Arrays.binarySearch(cumScoreDist, specProbThreshold);
//			if(index >= 0)
//				return cumScoreDist[index];
//			else
//			{
//				index = -index-1;
//				if(index > 0)
//					return cumScoreDist[index-1];
//				else
//					return 0;
//			}
			while(curIndex < cumScoreDist.length-1 && cumScoreDist[curIndex+1] <= specProbThreshold)
				++curIndex;
			
			return cumScoreDist[curIndex];
		}
		
		public double getSpecProb() {
			return specProb;
		}
		public double getPValue() {
			return pValue;
		}
		public String getResultStr() {
			return resultStr;
		}
		@Override
		public int compareTo(DBMatch arg0) {
			if(this.specProb < arg0.specProb)
				return -1;
			else if(this.specProb > arg0.specProb)
				return 1;
			else
				return 0;
		}
	}
}
