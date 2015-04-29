package Analysis;

import java.io.PrintStream;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.Spectrums.LabelledPeak;
import org.Spectrums.LinkedPeptide;
import org.Spectrums.Mass;
import org.Spectrums.Peak;
import org.Spectrums.PeakMassComparator;
import org.Spectrums.SimpleMatchingGraph;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumLib;
import org.Spectrums.SpectrumUtil;
import org.Spectrums.TheoreticalSpectrum;

/**
 * Linked peptide spectrum-spectrum matcher
 * Match Linked-peptide-spectrum with their unlinked counter part to evaluate
 * spectrum similarity
 * @author Jian
 *
 */
public class LPSSMatcher {
	double tolerance = 0.5;
	public LPSSMatcher(double tolerance){
		this.tolerance = tolerance;
	}
	
	
	public double projCosine(Spectrum linked, Spectrum unlinked, Spectrum unlinked2){
		TheoreticalSpectrum unlinkedTheo1 = new TheoreticalSpectrum(unlinked.peptide+unlinked.charge);
		TheoreticalSpectrum unlinkedTheo2 = new TheoreticalSpectrum(unlinked.peptide+unlinked.charge);
		LinkedPeptide linkedpep = new LinkedPeptide(linked.peptide, linked.charge);
		TheoreticalSpectrum linkedTheo = new TheoreticalSpectrum(linkedpep, (short)linked.charge, false);
		//projected-unlinked fragments
		SimpleMatchingGraph g1 = SpectrumUtil.constructMatchingGraph(unlinked, unlinkedTheo1, 0.5);
		SimpleMatchingGraph g2 = SpectrumUtil.constructMatchingGraph(unlinked2, unlinkedTheo2, 0.5);
		SimpleMatchingGraph g3 = SpectrumUtil.constructMatchingGraph(linked, linkedTheo, 0.5);
		
		
		
		return 0.0;
	}
	
	public Spectrum deCharge(Spectrum s){
		TheoreticalSpectrum th = null;
		int count=0;
		if(s.peptide.contains("--")){
			LinkedPeptide lp = new LinkedPeptide(s.peptide, s.charge);
			th = new TheoreticalSpectrum(lp, (short)s.charge, false);
		}else{
			th = new TheoreticalSpectrum(s.peptide+"."+s.charge);
		}
		Spectrum copy = new Spectrum(s);
		SimpleMatchingGraph g = SimpleMatchingGraph.getBipartiteMatching(copy, th, 0.5);
		Collection<Peak> newPeaks = new HashSet<Peak>();
		for(Iterator<Peak> it = copy.getPeak().iterator(); it.hasNext();){
			Peak curr = it.next();
			List<Peak> neighs  = g.getNeighbors(curr);
			if(neighs!=null && neighs.size() > 0){
				LabelledPeak lp = (LabelledPeak)neighs.iterator().next();
				if(lp.getCharge() > 1){
					Peak newPeak = new Peak(curr.getMass()*lp.getCharge()-Mass.PROTON_MASS*(lp.getCharge()-1), curr.getIntensity());
					it.remove();
					newPeaks.add(newPeak);
					count++;
				}
			}
		}
		//System.out.println("total decharged peaks: " + count);
		copy.getPeak().addAll(newPeaks);
		Collections.sort(copy.getPeak(), PeakMassComparator.comparator);
		return copy;
	}
	
	
	public static Map<String, List<Spectrum>> generateLPSMMap(){
		Map<String, List<Spectrum>> lpsmMap = new HashMap<String, List<Spectrum>>();
		return lpsmMap;
	}
	
	
	public static String[] getLinkPepSeq(String linkedPep){
			//System.out.println("linkedpep " + linkedPep);
			String[] tokens = linkedPep.split("\\||--");
			return new String[]{ Utils.StringUtils.getStrippedSeq(tokens[0]),
						Utils.StringUtils.getStrippedSeq(tokens[1])};		
	}
	
	public static void testRawProjLSSM(){
		String unlinkLib = "../mixture_linked/MXDB_SSMatch/unlinked.mgf";
		String linkLib = "../mixture_linked/MXDB_SSMatch/linked_fixed.mgf";
		SpectrumLib lib = new SpectrumLib(linkLib, "MGF");
		SpectrumLib lib2 = new SpectrumLib(unlinkLib, "MGF");
		toSeqKey(lib2);
		LPSSMatcher matcher = new LPSSMatcher(0.3);
		for(int i = 0; i < lib.getAllSpectrums().size(); i++){
			Spectrum linkSpect = (Spectrum)lib.getAllSpectrums().get(i);
			linkSpect.windowFilterPeaks2(10, 25);
			linkSpect.sqrtSpectrum();
			String[] pepseqs = getLinkPepSeq(linkSpect.peptide);
			//System.out.println(pepseqs[0] + "\t" + pepseqs[1]);
			for(int j = 0; j < pepseqs.length; j++){
				if(lib2.getSpectra(pepseqs[j]) != null){
					Spectrum s = lib2.getSpectra(pepseqs[j]).get(0);
					s.windowFilterPeaks2(10, 25);
					s.sqrtSpectrum();
					Spectrum decharge = matcher.deCharge(linkSpect);
					System.out.println(linkSpect.peptide + " & " + s.peptide + " projCos similarity: " +  s.projectedCosine(linkSpect, 0.333) + "\t" + linkSpect.projectedCosine(s, 0.333) 
							+ "\tafter-decharging\t" +  s.projectedCosine(decharge, 0.333) + "\t" + decharge.projectedCosine(s, 0.333));
				}
			}
		}
	}
	
	public static void toSeqKey(SpectrumLib lib){
		Map<String, List<Spectrum>> specMap = lib.getSpectrumLibrary();
		Map<String, List<Spectrum>> newMap = new HashMap<String, List<Spectrum>>();
		for(Iterator<String> it = specMap.keySet().iterator(); it.hasNext();){
			String pepKey = it.next();
			if(pepKey.matches(".*[\\+-0-9\\.]+.*")){
				String unModKey = pepKey.replaceAll("[0-9\\.\\+]", "");
				List<Spectrum> values = specMap.get(pepKey);
				newMap.put(unModKey, values);
				it.remove();
			}
		}
		specMap.putAll(newMap);
		lib.setSpectrumLibrary((Hashtable<String, List<Spectrum>>)specMap);
	}
	
	//fix the annotation in linked lib
	public static void fixLinkedLib(){
		String linkLib = "../mixture_linked/MXDB_SSMatch/linked.mgf";
		String fixLib = "../mixture_linked/MXDB_SSMatch/linked_fixed.mgf";
		SpectrumLib lib = new SpectrumLib(linkLib, "MGF");
		PrintStream out = Utils.FileIOUtils.getOutStream(fixLib);
		for(int i = 0; i < lib.getAllSpectrums().size(); i++){
			Spectrum s = (Spectrum)lib.getAllSpectrums().get(i);
			String linkseq=s.spectrumName.split("PROTEIN: ")[1].split(" ")[0];
			LinkedPeptide lp = new LinkedPeptide(linkseq, 1);
			System.out.println(linkseq + "\tparsed:\t" + lp);
			int charge = LinkedPeptide.guessCharge(lp, s.parentMass, 3, 6);
			System.out.println(lp.getParentmass() + "\t" + s.parentMass + "\tcharged:\t" + charge);
			s.peptide = lp.toString();
			s.charge = charge;
			out.println(s);
		}
		
	}
	
	public static void main(String[] args){
		//fixLinkedLib();
		testRawProjLSSM();
	}
	
	
	
	
}
