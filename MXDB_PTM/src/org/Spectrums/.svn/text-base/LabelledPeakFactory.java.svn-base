package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

public class LabelledPeakFactory {
	private static int currentIndex = 0;
	private static List<LabelledPeak> peakList = new ArrayList<LabelledPeak>();
	private static boolean poolMode = false; //determine we should recycle the objects
	                                         //set to true if expect a lot of temporary creation of peaks
	public static boolean isPoolMode() {
		return poolMode;
	}

	public static void setPoolMode(boolean poolMode) {
		LabelledPeakFactory.poolMode = poolMode;
	}

	public static LabelledPeak createLabelledPeak(double moz, double intensity, String type, short pos, short charge){
		if(!poolMode){
			LabelledPeak lp = new LabelledPeak(moz, intensity, type, pos, charge);
			return lp;
		}
		if(currentIndex < peakList.size()){
			LabelledPeak lp = peakList.get(currentIndex);
			lp.setMoz(moz);
			lp.setIntensity(intensity);
			lp.setType(type);
			lp.setPos(pos);
			lp.setCharge(charge);
			currentIndex++;
			return lp;
		}else{
			LabelledPeak lp = new LabelledPeak(moz, intensity, type, pos, charge);
			peakList.add(lp);
			return lp;
		}
		
	}
	
	public static void resetFactory(){
		currentIndex = 0;
	}
	
	public static void clearFactory(){
		peakList.clear();
	}
}
