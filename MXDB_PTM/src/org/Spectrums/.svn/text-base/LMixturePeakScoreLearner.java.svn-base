package org.Spectrums;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

public class LMixturePeakScoreLearner implements PeakComparator, Serializable{
	private static final long serialVersionUID = 193712381028301L;
	private MixturePeakScoreLearner core;
	public LMixturePeakScoreLearner(String trainingFile){
		this.core = new MixturePeakScoreLearner(trainingFile);
		this.core.getMixtureIonCount();
	}
	@Override
	public double compare(Peak p1, Peak p2) {
		
		MixturePeak lp = (MixturePeak)p1;
		if(p1 == null){
			return 0.0;
		}else{
			int pepCharge = lp.getPep().getCharge();
			int oldParentCharge = lp.getParent().charge;
			if(lp.getParent().charge == 2){
				lp.getPep().setCharge((short)2);
			}
			if(lp.getParent().charge <= 3){
				lp.getPep().setCharge((short)2);
			}else{
				lp.getPep().setCharge((short)3);
			}
			int newParentCharge = oldParentCharge > 6 ? 6 : oldParentCharge;
			newParentCharge = oldParentCharge <= 3 ? oldParentCharge+1 : oldParentCharge;
			lp.getParent().charge = newParentCharge;
			
			double score = 0.0;
			if(TheoreticalSpectrum.isLinkedPeak(lp.getPep(), lp)){
				int oldCharge = lp.getCharge();
				int newCharge = LinkedPeptide.transformPeakCharge(lp.getCharge(), lp.getParent().charge);
				if(newCharge > 0){
					lp.setCharge((short)newCharge);
					//System.out.println("comparing " + lp);
					score = this.core.compare(p1, p2);
					lp.setCharge((short)oldCharge);
				}
			}else{
				if(lp.getCharge() <= lp.getPep().getCharge()){
					score = this.core.compare(p1, p2);
					
				}
			}
			lp.getParent().charge = oldParentCharge;
			lp.getPep().setCharge((short)pepCharge);
			return score;
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
	public static LMixturePeakScoreLearner loadComparator(String file){
		try{
			BufferedInputStream bi = new BufferedInputStream(new FileInputStream(file));
			ObjectInputStream oi = new ObjectInputStream(bi);
		    Object o = oi.readObject();
		    return (LMixturePeakScoreLearner)o;
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	
}
