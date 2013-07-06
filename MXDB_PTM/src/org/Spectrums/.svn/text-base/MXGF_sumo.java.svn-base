package org.Spectrums;

public class MXGF_sumo extends MXGF{
	double resolution=1.0;
	private double[] prefixScores;
	private double[] suffixScores;
	private double[] tagScores;
	private double[][] prefixProb;
	private double[][] suffixProb;
	private double[][] tagProb;
	private double maxPrefixScore=this.maxMass;
	private double maxSuffixScore=this.maxMass;
	private double maxTagScore=this.maxMass;
	private double minPrefixScore=this.minScore;
	private double minSuffixScore=this.minScore;
	private double minTagScore=this.minScore;
	private int[] validSiteAA; //linking site at K
	private int[] validTagAA;
	private int[] substrateMasses;
	private int[] tagMasses;

	public MXGF_sumo(double precursor, double[] prefix, double[] suffix, double[] tags){
		this.maxMass = precursor; //allow some extra room for max mass
		this.setScore(prefix, suffix, tags);
	}
	
	public void setScore(double[] preffix, double[] suffix, double[] tags){
		this.prefixScores = preffix;
		this.suffixScores = suffix;
		this.tagScores = tags;
	}
	
	public void setMasses(int[] substrateMasses, int[] tagMasses){
		this.substrateMasses = substrateMasses;
		this.tagMasses = tagMasses;
	}
	
	public void setLinkerMasses(int[] linkerMasses){
		this.validSiteAA = linkerMasses;
	}
	
	public void initializeSUMO(){
		double maxScore, minScore;
		this.validSiteAA = new int[]{128+599}; //linking site at K
		this.validTagAA = new int[this.validAAMass.length+1];
		//allowing pyro-Q for tag
		for(int i = 0; i < this.validAAMass.length; i++){
			this.validTagAA[i] = this.validAAMass[i];
		}
		this.validTagAA[this.validTagAA.length-1]=111;
		for(int i = 0; i < 700; i++){
			//this.suffixScores[i] = 0;
		}

		int length = (int)Math.ceil(maxMass/resolution)+1;
		this.maxPrefixScore = this.computeMaxScore(prefixScores);
		this.minPrefixScore = this.computeMinScore(prefixScores);
		this.prefixProb = new double[length][(int)Math.ceil(this.maxPrefixScore-this.minPrefixScore+1)];
		
		this.maxSuffixScore = this.computeMaxScore(suffixScores) + this.maxPrefixScore;
		this.minSuffixScore = this.computeMinScore(suffixScores) + this.minPrefixScore;
		this.suffixProb = new double[length][(int)Math.ceil(this.maxSuffixScore-this.minSuffixScore+1)];
		
		this.maxTagScore = this.computeMaxScore(tagScores);
		this.minTagScore = this.computeMinScore(tagScores);
		this.tagProb = new double[length][(int)Math.ceil(this.maxTagScore - this.minTagScore + 1)];
		
	}
	
	public void computeProb(){
		
		//masking tag prm position in substrate score vectors
		for(int i = 0; i < this.tagMasses.length; i++){
			prefixScores[tagMasses[i]] = 0;
			suffixScores[tagMasses[i]] = 0;
		}
		
		//computing prefix prob
		prefixProb[0][(int)(-1*this.minPrefixScore)] = 1;
		int minS = (int)Math.floor(this.minPrefixScore);
		//System.out.println("offset: " + scoreOffSet);
		System.out.println("max score bin: " + this.prefixScores.length);
		System.out.println("max mass bin: " + this.maxMass);
		for(int m1 = 0; m1 < maxMass; m1++){
			for(int s = minS; s < this.maxPrefixScore; s++){
				for(int a = 0; a < this.validAAMass.length; a++){
					int prevMass = m1 - this.validAAMass[a];
					int prevScore = (int)(s - this.prefixScores[m1]);
 					double prevCount = 0;
					if(prevMass >= 0 && prevScore > this.minPrefixScore && prevScore < this.maxPrefixScore){
						prevScore = prevScore - minS;
						prevCount = this.prefixProb[prevMass][prevScore];
						this.prefixProb[m1][s - minS] += prevCount*0.05;//*this.aaFreq[a];
					}
				}
			}
		}
		//computing suffix prob
		//suffixProb[0][(int)(-1*this.minSuffixScore)] = 1;
		minS = (int)Math.floor(minSuffixScore);
		suffixProb[this.validSiteAA[0]][(int)(this.tagScores[this.validSiteAA[0]])-minS] = 1;
		int minS2 = (int)Math.floor(minPrefixScore);
		for(int m1 = this.validSiteAA[0]; m1 < maxMass; m1++){
			for(int s = minS; s < this.maxSuffixScore; s++){
				//loop through all aa in suffix
				for(int a = 0; a < this.validAAMass.length; a++){
					int prevMass = m1 - this.validAAMass[a];
					int prevScore = (int)(s - this.suffixScores[m1]);
 					double prevCount = 0;
					if(prevMass >= 0 && prevScore > this.minSuffixScore && prevScore < this.maxSuffixScore){
						prevScore = prevScore - minS;
						prevCount = this.suffixProb[prevMass][prevScore];
						this.suffixProb[m1][s - minS] += prevCount*0.05;//*this.aaFreq[a];
					}
				}
				//jumping from prefix
				for(int a = 0; a < this.validSiteAA.length; a++){
					//System.out.println("jumping mass is: " + this.validSiteAA[a]);
					int prevMass = m1 - this.validSiteAA[a];
					int prevScore = (int)(s - this.suffixScores[m1]);
 					double prevCount = 0;
					if(prevMass >= 0 && prevScore > this.minPrefixScore && prevScore < this.maxPrefixScore){
						//System.out.println("jumping from prefix: " + a);
						prevScore = prevScore - minS2;
						prevCount = this.prefixProb[prevMass][prevScore];
						this.suffixProb[m1][s - minS] += prevCount*0.05;//*this.aaFreq[a];
					}
				}
			}
			
			
		}
		
		//computing tag prob
		tagProb[0][(int)(-1*this.minTagScore)] = 1;
		minS = (int)Math.floor(minTagScore);
		//System.out.println("offset: " + scoreOffSet);
		for(int m1 = 0; m1 < maxMass; m1++){
			for(int s = minS; s < this.maxTagScore; s++){
				for(int a = 0; a < this.validTagAA.length; a++){
					int prevMass = m1 - this.validTagAA[a];
					int prevScore = (int)(s - this.tagScores[m1]);
 					double prevCount = 0;
					if(prevMass >= 0 && prevScore > minTagScore && prevScore < maxTagScore){
						prevScore = prevScore - minS;
						prevCount = this.tagProb[prevMass][prevScore];
						this.tagProb[m1][s - minS] += prevCount*0.05;//*this.aaFreq[a];
					}
				}
			}
		}
		
	}
	
	public double[] getSUMOProb(double substrateScore, double tagScore, double mass){
		double substrateP = getSpecProb(substrateScore, this.maxSuffixScore, this.minSuffixScore, mass, this.suffixProb);
		double tagP = getSpecProb(tagScore, this.maxTagScore, this.minTagScore, mass, this.tagProb);
		System.out.println("current tagscore: " + tagScore + "\t" + "maxTagScore\t" + this.maxTagScore);
		System.out.println("current mass: " + mass + "\t" + "maxmass\t" + this.maxMass);
		return new double[]{substrateP, tagP};
	}
	
	private double getSpecProb(double score, double maxScore, double minScore, double mass, double[][] probTable){
		int s = (int)Math.round(score-minScore);
		int m1 = (int)Math.round(this.scaleMass(mass));
		//System.out.println("m1 " + m1 + "\t" + s);
		double total = 0;
		double higher = 0.0;
		double highest = 0.0;
		for(int currentS = 0; currentS < maxScore-minScore; currentS++){
			if(currentS >= s){
				higher+=probTable[m1][currentS];
			}
			total += probTable[m1][currentS];
			if(probTable[m1][currentS] > 0){
				highest=probTable[m1][currentS];
			}
			//System.out.println((currentS + minScore) + "\t" + table[m1][index][currentS]);
		}
		if(higher == 0){
			higher=highest; //if higher is zero we use highest
		}
		//System.out.println("total possible candidates: " + this.pepCounts[m1]);
		//System.out.println("total: " + total + "\t" + "higher " + higher);
		return higher;///total;
	}
	
	
	public static void testSUMOProb(){	
		
	}
	
	
	
	
}
