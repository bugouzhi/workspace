package org.Spectrums;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.collections.MultiMap;
import org.apache.commons.collections.map.MultiValueMap;
/**
 * Peptidic features in LC/MS maps
 * @author Jian Wang
 *
 */
public class MSFeature {
	private String id;
	private double quality;
	private double intensity;
	private double mz;  
	private double rt; 
	private MultiValueMap peakList;
	private int charge;
	private double minRT;
	private double maxRT;
	private int Scan;
	public int getScan() {
		return Scan;
	}
	public void setScan(int scan) {
		Scan = scan;
	}
	public int getCharge() {
		return charge;
	}
	public void setCharge(int charge) {
		this.charge = charge;
	}
	public MSFeature(){
		this.id = "";
		this.quality = 0.0;
		this.intensity = 0.0;
		this.mz = 0.0;
		this.rt = 0.0;
		this.peakList = null;
	}
	
	public boolean isWithinFeature(double monoMass, double RTtime, int charge, double massTolerance){
		return (RTtime >= this.minRT && RTtime <= this.maxRT && charge==this.charge && Math.abs(monoMass-this.mz) < massTolerance);			
	}
	
	public boolean isWithnFeature(double precursorMass, double RTtime, double massTolerance){
		if(!(RTtime >= this.minRT && RTtime <= this.maxRT)){
			return false;
		}
		for(Iterator it = this.peakList.values().iterator(); it.hasNext();){
			double mass = (Double)it.next();
			//System.out.println("feature mass: " + mass + "\t" + "ms2 mass: " + precursorMass);
			if(Math.abs(mass-precursorMass) < massTolerance){
				return true;
			}
		}
		return false;
	}
	
	public boolean isMatchFeature(Peptide p, double massTolerance){
		double shiftmass = this.mz - (this.mz*100/1000000);
		return(p.getCharge() == this.getCharge()
				&& Math.abs(p.getParentmass() - shiftmass ) < massTolerance);
//		return (Math.abs(p.getParentmass() - this.mz) < massTolerance);
	}

	public boolean isMatchFeature(Peptide p, double massTolerance, double RTtime){
		double shiftmass = this.mz - (this.mz*100/1000000);
		return( RTtime >= this.minRT && RTtime <= this.maxRT
				&& p.getCharge() == this.getCharge()
				&& Math.abs(p.getParentmass() - shiftmass ) < massTolerance);
//		return (Math.abs(p.getParentmass() - this.mz) < massTolerance);
	}
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public double getQuality() {
		return quality;
	}
	public void setQuality(double quality) {
		this.quality = quality;
	}
	public double getIntensity() {
		return intensity;
	}
	public void setIntensity(double intensity) {
		this.intensity = intensity;
	}
	public double getMz() {
		return mz;
	}
	public void setMz(double mz) {
		this.mz = mz;
	}
	public double getRt() {
		return rt;
	}
	public void setRt(double rt) {
		this.rt = rt;
	}
	public MultiValueMap getFeatureList() {
		return peakList;
	}
	public void setFeatureList(MultiValueMap featureList) {
		this.peakList = featureList;
		List<Double> keyList = new ArrayList<Double>();
		keyList.addAll(featureList.keySet());
		Collections.sort(keyList);
		this.minRT = keyList.get(0);
		this.maxRT = keyList.get(keyList.size()-1);
	}
	
//	public String toString(){
//		StringBuffer buff = new StringBuffer();
//		buff.append("id: "+ "\t" + this.id + "\n");
//		buff.append("quality: "+ "\t" + this.quality + "\n");
//		buff.append("charge: "+ "\t" + this.charge + "\n");
//		buff.append("intensity: "+ "\t" + this.intensity + "\n");
//		buff.append("rep retention: " + "\t" + this.rt + "\n");
//		buff.append("rep mz: " + "\t" + this.mz + "\n");
//		List<Double> rtList = new ArrayList(); 
//		if(this.peakList != null){
//			rtList.addAll(peakList.keySet());
//			Collections.sort(rtList);
//			for(int i = 0; i < rtList.size(); i++){
//				buff.append("RT:\t" + rtList.get(i) + "\tmasses:\t");
//				Iterator it = this.peakList.getCollection(rtList.get(i)).iterator();
//				while(it.hasNext()){
//					buff.append(it.next() + "\t");
//				}
//				buff.append("\n");
//		
//			}
//		}
//		return buff.toString();
//	}

	public String toString(){
		StringBuffer buff = new StringBuffer();
		buff.append("feature: " +this.id + "\t");
		buff.append(this.quality + "\t");
		buff.append(this.charge + "\t");
		buff.append(this.intensity + "\t");
		buff.append(this.rt + "\t");
		buff.append(this.Scan + "\t");
		buff.append(this.mz + "\t");
//		List<Double> rtList = new ArrayList(); 
//		if(this.peakList != null){
//			rtList.addAll(peakList.keySet());
//			Collections.sort(rtList);
//			for(int i = 0; i < rtList.size(); i++){
//				buff.append("RT:\t" + rtList.get(i) + "\tmasses:\t");
//				Iterator it = this.peakList.getCollection(rtList.get(i)).iterator();
//				while(it.hasNext()){
//					buff.append(it.next() + "\t");
//				}
//				//buff.append("\n");
//		
//			}
//		}
		return buff.toString();
	}

	

}
