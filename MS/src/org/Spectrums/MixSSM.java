package org.Spectrums;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import IO.MZXMLReader;

/**
 * Describe a mixture or multiplexed spectrum-spectrum matches
 * @author Jian
 *
 */
public class MixSSM {
	private Spectrum query;
	private Collection<Spectrum> matches;
	public MixSSM(Spectrum query, Collection<Spectrum> matches){
		this.query = query;
		this.matches = matches;
	}
	public Spectrum getQuery() {
		return query;
	}
	public void setQuery(Spectrum query) {
		this.query = query;
	}
	public Collection<Spectrum> getMatches() {
		return matches;
	}
	public void setMatches(Collection<Spectrum> matches) {
		this.matches = matches;
	}
	
	public List<Peak> getShareMatchedPeaks(Spectrum match){
		return getShareMatchSpect(match).getPeak();
	}
	
	public Spectrum getShareMatchSpect(Spectrum match){
		List<Peak> unique = new ArrayList<Peak>();
		List<Peak> matchedPeaks = new ArrayList<Peak>();
		Spectrum otherMatches = new Spectrum();
		for(Iterator<Spectrum> it = matches.iterator(); it.hasNext();){
			Spectrum curr = it.next();
			if(curr != match){
				SpectrumUtil.mergeInto(otherMatches, curr, 0.05);
			}
		}
		Collections.sort(otherMatches.getPeak(), PeakMassComparator.comparator);
		Spectrum proj = query.project(otherMatches, 0.05);
		//System.out.println("proj-size: " + proj.getPeak().size());
		return proj;	
	}
	
	public Spectrum getUniqueSpect(Spectrum match){
		List<Peak> unique = new ArrayList<Peak>();
		List<Peak> matchedPeaks = new ArrayList<Peak>();
		Spectrum otherMatches = new Spectrum();
		for(Iterator<Spectrum> it = matches.iterator(); it.hasNext();){
			Spectrum curr = it.next();
			if(curr != match){
				SpectrumUtil.mergeInto(otherMatches, curr, 0.05);
			}
		}
		Collections.sort(otherMatches.getPeak(), PeakMassComparator.comparator);
		Spectrum proj= match.removeSharePeaks(otherMatches, 0.05);
		//System.out.println("proj-size: " + proj.getPeak().size());
		return proj;	
	}
	
	
	public static MixSSM createMixSSM(int queryScan, MZXMLReader reader, SpectrumLib lib, List<String> matches){
		Spectrum query = reader.getSpectrum(queryScan);
		return createMixSSM(query, lib, matches);
	}
	
	public static MixSSM createMixSSM(Spectrum query, SpectrumLib lib, List<String> matches){
		if(matches == null){
			return null;
		}
		List<Spectrum> matchSpects = new ArrayList<Spectrum>(matches.size());
		for(int j = 0; j < matches.size(); j++){
			String pepKey1 = matches.get(j);
			Spectrum libSpect = lib.getSpectra(pepKey1).get(0);
			libSpect.mergePeaks(libSpect, 0.05);
			libSpect.filterPeaks(30);
			matchSpects.add(libSpect);
		}
		return new MixSSM(query, matchSpects);
	}
	
}
