package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

/**
 * A MSFeature with acquisition information
 * @author Jian
 *
 */
public class MSFeatureAcq extends MSFeature{
	private List<Spectrum> acqInfo;
	
	public MSFeatureAcq(){
		super();
		this.acqInfo = new ArrayList<Spectrum>();
	}
	
	public MSFeatureAcq(MSFeature feature){
		this.setId(feature.getId());
		this.setQuality(feature.getQuality());
		this.setMz(feature.getMz());
		this.setIntensity(feature.getIntensity());
		this.setMinRT(feature.getMinRT());
		this.setMaxRT(feature.getMaxRT());
		this.acqInfo = new ArrayList<Spectrum>();
	}
	
	public List<Spectrum> getAcqInfo() {
		return acqInfo;
	}

	public void setAcqInfo(List<Spectrum> acqInfo) {
		this.acqInfo = acqInfo;
	}
	
}
