package org.Spectrums;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeModel;

import org.systemsbiology.jrap.stax.*;
import Utils.FileIOUtils;
//a combinatoric peptide library, represent by some constant
//positions and some variable positions

public class CombinatoryPeptides {
	private String combPeptide; //the master peptide pattern, use to generate the combinatoric library
	private List<String> peptides;
	public CombinatoryPeptides(String peptidePattern){
		this.combPeptide = peptidePattern;	
	}
	
	public List<String> generateAllPeptides(){
		this.peptides = new ArrayList<String>();
		StringBuffer buff = new StringBuffer();
		List<StringBuffer> current = new ArrayList<StringBuffer>();
		current.add(buff);
		List<StringBuffer> peptides = generatePeptide(0, current);
		List<String> peps = new ArrayList<String>();
		for(int i = 0; i < peptides.size(); i++){
			StringBuffer pep = peptides.get(i);
			double pm = PeptideMassAnalysis.computeMolecularMass(pep.toString());
			//pep.insert(1, "+42.010564686");
			peps.add(peptides.get(i).toString());
			System.out.println(peptides.get(i) + "\tmass: " + pm);
		}
		System.out.println("Total peptides generated: " + peptides.size());
		return peps;
	}
	
	public List<String> generateDigestedPeptides(){
		this.peptides = new ArrayList<String>();
		StringBuffer buff = new StringBuffer();
		List<StringBuffer> current = new ArrayList<StringBuffer>();
		current.add(buff);
		List<StringBuffer> peptides = generatePeptide(0, current);
		Set<String> digested = new HashSet<String>();
		for(int i = 0; i < peptides.size(); i++){
			digested.addAll(this.cutPeptide(peptides.get(i).toString()));
		}
		System.out.println("Total peptides generated: " + digested.size());
		List<String> digestedSet = new ArrayList<String>();
		digestedSet.addAll(digested);
		return digestedSet;
	}
	
	public List<String> generateAllLinkedPeptides(){
		List<String> peptides, linkedPeptides;
		linkedPeptides = new ArrayList<String>();
		peptides = this.generateAllPeptides();
		for(Iterator<String> it = peptides.iterator(); it.hasNext();){
			String current = it.next();
			linkedPeptides.addAll(crossLinkPeptide(current));
			int index = linkedPeptides.size()-2;
			Peptide linked = new LinkedPeptide(linkedPeptides.get(index),1);
			Peptide linked2 = new LinkedPeptide(linkedPeptides.get(index+1),1);
			double p = PeptideMassAnalysis.computeMbyZ(current, 1);
			System.out.println("peptide: " + current +  " mass: " + (p + 42.010564686) + "\t"
					+ " linked: " + linkedPeptides.get(index) 
					+ " mass: " + linked.getParentmass() + "\t"
					+ linkedPeptides.get(index+1) + " mass: " + linked2.getParentmass());
		}
		return linkedPeptides;
	}
	
	public List<String> generateAllLinkedPairs(){
		List<String> linkedPeptides = new ArrayList<String>();
		List<String> peptides = this.generateAllPeptides();
		peptides = this.addMod(peptides, 42.010564686);
		List<String> digested = this.cutPeptides(peptides);
		System.out.println("We have " + digested.size() + " unique digested peptides");
		linkedPeptides = this.crossLinkPeptides(digested);
		System.out.println("We cross-link " + linkedPeptides.size() + " unique linked peptides");
		return linkedPeptides;
	}
	
	public List<StringBuffer> generatePeptide(int i, List<StringBuffer> current){
		if(i == this.combPeptide.length()-1){
			//System.out.println("0");
			return append(current, combPeptide.charAt(i));
		}else if(this.combPeptide.charAt(i) == '['){
			//System.out.println("1");
			List<StringBuffer> newSet = new ArrayList<StringBuffer>();
			int j = i+1;	
			for(; j < this.combPeptide.length() 
				&& this.combPeptide.charAt(j) != ']'; j++){
				newSet.addAll(append(current, combPeptide.charAt(j)));
			}
			if(j == this.combPeptide.length()-1){
				return newSet;
			}else{
				return generatePeptide(j+1, newSet);
			}
		}else{
			//System.out.println("2");
			List<StringBuffer> newSet = new ArrayList<StringBuffer>();
			newSet.addAll(append(current, combPeptide.charAt(i)));
			return generatePeptide(i+1, newSet);
		}
	}
	
	//cut the peptide into cross-llnked pairs
	public List<String> crossLinkPeptide(String peptide){
		List<String> peptides = new ArrayList<String>();
		peptides.add(peptide.substring(2, 10) + "--" + peptide.substring(14));
		peptides.add(peptide.substring(2, 14) + "--" + peptide.substring(14));
		return peptides;
	}
	
	public List<String> crossLinkPeptides(List<String> peptides){
		List<String> crosslinks = new ArrayList<String>();
		Set<String> crosslinkSet = new HashSet<String>();
		for(int i = 0; i < peptides.size(); i++){
			for(int j=i+1; j < peptides.size(); j++){
				String crosslink = peptides.get(j) + "--" + peptides.get(i);
				if(!crosslinkSet.contains(crosslink)){
					crosslinks.add(peptides.get(j) + "--" + peptides.get(i));
					System.out.println("cross-linked peptide: " + crosslinks.get(crosslinks.size()-1));
				}
			}
		}
		crosslinks.addAll(crosslinkSet);
		return crosslinks;
	}
	
	public List<String> cutPeptide(String peptide){
		List<String> peptides = new ArrayList<String>();
		peptides.add(peptide.substring(0, 10));
		peptides.add(peptide.substring(14));
		peptides.add(peptide.substring(0, 14));
		peptides.add(peptide.substring(14));
		return peptides;
	}
	
	public List<String> cutPeptides(List<String> peptides){
		List<String> digested = new ArrayList<String>();
		for(int i = 0; i < peptides.size(); i++){
			int ind = -1;
			String pep = peptides.get(i);
			while(ind < pep.length()-1){
				digested.add(pep.substring(ind+1));
				ind = pep.indexOf('R', ind+1);
			}
		}
		return digested;
	}
	
	private List<StringBuffer> append(List<StringBuffer> current, char c){
		List<StringBuffer> newStr = new ArrayList<StringBuffer>();
		for(Iterator<StringBuffer> it = current.iterator(); it.hasNext();){
			StringBuffer buffer = new StringBuffer();
			buffer.append(it.next());
			buffer.append(c);
			newStr.add(buffer);
		}
		return newStr;
	}
	
	//add modification to generated peptides
	public static List<String> addMod(List<String> peps, double mod){
		List<String> modified = new ArrayList<String>();
		for(int i = 0; i < peps.size(); i++){
			StringBuffer buff = new StringBuffer(peps.get(i));
			buff.insert(1, "+"+mod);
			modified.add(buff.toString());
		}
		return modified;
	}
	
	//same as addMod by keep the original peptides
	public static List<String> appendMod(List<String> peps, double mod){
		List<String> modified = new ArrayList<String>();
		for(int i = 0; i < peps.size(); i++){
			StringBuffer buff = new StringBuffer(peps.get(i));
			buff.insert(1, "+"+mod);
			modified.add(peps.get(i));
			modified.add(buff.toString());
		}
		return modified;
	}
	
	public static List<String> addSubSequence(List<String> peptides, int maxBegin, int minEnd){
		List<String> subPeps = new ArrayList<String>();
		for(int i = 0; i < peptides.size(); i++){
			String currentPep = peptides.get(i);
			for(int begin = 0; begin < maxBegin; begin++){
				for(int end = 0; end < minEnd; end++){
					String subPep = currentPep.substring(begin, end);
					subPeps.add(subPep);
				}
			}
		}
		return subPeps;
	}
	
	public static void testGeneratePeptide(){
		String pattern = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCC[EA]KQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYF[YL]APELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVL[AT]SSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKEC[CK]D[KP][PC]LLEKSHCIAEVEKDAIPE[ND]LPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACY[TS][ST]VFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQN[AE]LIVRYTR[KR]VPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPC[TA]EDYLSLILNRLCVLHEKTPV[SE][ES]KVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFV[DG]KCCAADDKEACFAVEGPKLVVSTQTALA";
		CombinatoryPeptides combPeps = new CombinatoryPeptides(pattern);
		combPeps.generateAllPeptides();
		//combPeps.generateAllLinkedPeptides();
	}
	
	public static void testMatchPeptide(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("TR[TN][PG][AY]K[EF][IQ][DS]RASSR[TN][PG][AY]K[EF][IQ][DS]R");
		List<String> peptides = combPeps.generateAllPeptides();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(6);
		factory.indexPeptideByParentMass(0.01);
		String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\20100215_nLCFTMSMS_CLPL_2pmol_1841_MSonly.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		System.out.println("total library spectrum: " + specList.size());
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.scanNumber > 3300 || s.parentMass < 450){
				continue;
			}else{
				System.out.println(s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
			}
			//s.windowFilterPeaks(6, 25);
//			s.computePeakRank();
			SpectrumLib lib  = factory.createCandidateSpectrumLibX(s, 0.01, false);
			List<Spectrum> candidates = lib.getSpectrumList();
			System.out.println(s.spectrumName + " has " + lib.getAllSpectrums().size() + " candidates");
//			for(int j = 0; j < candidates.size(); j++){
//				TheoreticalSpectrum t = (TheoreticalSpectrum)candidates.get(j);
//				double explained = t.explainedPeaks2(t, s, 0.1);
//				System.out.println(s.spectrumName + " match to candidate: " + candidates.get(j).peptide + " explained intensity: " + explained);
//				
//			}
			//SpectrumLibSearcher searcher = new SpectrumLibSearcher(lib.getAllSpectrums(),  scorer, scorer);
			//searcher.topSpectra(s, 3);
		}
	}
	
	public static CandidateSpectrumLibFactory indexCrossLinkedPeptides(List<String> peptides, int minCharge, int maxCharge){
		//List<String> peptides = Utils.FileIOUtils.createListFromFile("..\\mixture_linked\\library_14mer_albumin_peptides.txt");
		//peptides = addMod(peptides, 42.010564686);
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);		
		factory.setMinCharge(1);
		factory.setMaxCharge(1);
		factory.indexPeptideByParentMass(0.5);
		factory.crossLinkAllPeptides(minCharge, maxCharge);	
		return factory;
	}
	
	public static void testMatchPeptideMasses(){
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("TR[TN][PG][AY]K[EF][IQ][DS]R");
		//CombinatoryPeptides combPeps2 = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		CombinatoryPeptides combPeps = new CombinatoryPeptides("TR[DW][GH]TGCA[EA]VKL[FT]D[DY]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[DS][IQ][EF]K[AY][PG][TN]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[YG]A[EQ]V[EL]LKT[AL]G[FD]V[YA]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]LKT[AL]G[FD]V[YA]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]VKL[AL]G[FD]T[YA]R"); 
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("KA[ARDEHLKMFPSTYV]D[ARDEHLKMFPSTYV]ES[ARDEHLKMFPSTYV]LRAK"); 
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("TALH[ARDEHLKMFPSTYV]K[ARDEHLKMFPSTYV]S[ARDEHLKMFPSTYV]TFR");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("A[FILVYW]K[ARDEHLKMFPSTYV][DE]T[ARDEHLKMFPSTYV]FRAK"); 
		///List<String> peptides = combPeps.generateAllPeptides();
		//peptides.addAll(combPeps2.generateAllPeptides());
		//List<String> peptides = new ArrayList<String>();
		List<String> peptides = Utils.FileIOUtils.createListFromFile("..\\mixture_linked\\database\\lib_ucsd3_plus_decoy.txt");
		//peptides.addAll(ecoPep);
		//CombinatoryPeptides combPeps2 = new CombinatoryPeptides("QQQTGG"); 
		//List<String> peptides2 = combPeps2.generateAllPeptides();
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);	
		factory.indexPeptideByParentMass(0.01);
		factory.insertPTM(42.010565, 1, 1);
//		factory.insertPTM(-57, new char[]{'C'}, 1);
//		CandidateSpectrumLibFactory factory2 = 	
//			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides2);
//		factory2.setMinCharge(1);
//		factory2.setMaxCharge(1);	
//		factory2.indexPeptideByParentMass(0.1);
		CrossLinker CROSSLINKER = new CrossLinker(Mass.DSSLINKER_MASS, new int[]{CrossLinker.ANYPOSITION}, new int[]{CrossLinker.ANYPOSITION}, new char[]{'C'}, new char[]{'C'});
		//factory.crossLinkAllPeptides(factory, 4, 6, 7, 7);
		factory.crossLinkAllPeptides(factory, 3, 5, CROSSLINKER);
		//TheoreticalSpectrum.prefixIons = Mass.standardPrefixesX;
		//TheoreticalSpectrum.suffixIons = Mass.standardSuffixesX;
		
		Mass.DSSLINKER_MASS = -116.0430; //138.06;//(Mass.C13-Mass.C12)*6;
		factory.setMatchCharge(true);
		String training2 = "..\\MSPLib\\Lib\\yeast.msp";
		//String training2 = "..\\mixture_linked\\mixtures100000_alpha0.5.mgf";
        //SpectrumComparator comp2 = SpectrumUtil.getLPeakRankBaseScorer(training2);
		SpectrumComparator comp2 = null;
        //SpectrumComparator comp2 = SpectrumUtil.getRankBaseScorer(training2);
		String spectrumLibFile = "..\\mixture_linked/linked_peptide_library/toni/110617_Crosslink/tk110610_Nuno_DMSO_peptides_6.mzXML";
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\/lib_ucsd/lib_ucsd3/11_0415_pep3_ox_x5.mzxml";		
//		PrecursorMassChecker checker = new PrecursorMassChecker(spectrumLibFile);
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		//List<Spectrum> specList = reader.readAllMS2Spectra();
		//reader = null;	
//		List<Spectrum>specList = lib.getAllSpectrums();
//		CandidateSpectrumLibFactory factory2 = 
//			CandidateSpectrumLibFactory.createFactoryFromPeptide("..\\mixture_linked\\linked_peptides.txt");
		PreciseCandidatesFactory pFactory = new PreciseCandidatesFactory(factory, spectrumLibFile);
		for(;reader.hasNext();){
			Spectrum s = reader.next();
			System.out.println("number of peaks: " + s.getPeak().size());
			s.windowFilterPeaks(8, 25);
			s.computePeakRank();
			if(s.scanNumber != 2718){
				//continue;
			}
//			if(s.scanNumber > 136000 || s.parentMass < 350 ){
	
//			}else{
				System.out.println("Query: " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
//			}
//			List<Peptide> preciseCandidates = pFactory.getCandidateByMass(s, 30);
			List<Peptide> preciseCandidates = pFactory.getCandidateByMassCrude(s, 2);
			for(int j = 0; j < preciseCandidates.size(); j++){
				Peptide p = preciseCandidates.get(j);
				//TheoreticalSpectrum t = generateLinkedTheoreticalSpectrum((LinkedPeptide)p);
				TheoreticalSpectrum t = new TheoreticalSpectrum((LinkedPeptide)p, p.getCharge(), false);	
				matchMS2(t, s, comp2, 0.5);
			}
		}
	}
	
	
	public static void testMatchSingle(){
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]LKT[AL]G[FD]V[YA]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("KEYPTENYTDDPDC");
		CombinatoryPeptides combPeps = new CombinatoryPeptides("CDPDDTYNETPYEK");
		List<String> peptides = combPeps.generateAllPeptides();
		//peptides = addMod(peptides, 42.010564686);
		//peptides = combPeps.cutPeptides(peptides);
		System.out.println("we have peptides: " + peptides.size());
		//peptides = addMod(peptides);
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(4);
		factory.indexPeptideByParentMass(0.1);
		//factory.reIndexPeptideByParentMass(0.1);
		//factory.insertPTM(Mass.DSSDANGLE_MASS, new char[]{'K', 'Y', 'Q'}, 1);
		factory.insertPTM(1826.7, new char[]{'K'}, 1);
		factory.insertPTM(-57, new char[]{'C'}, 1);
		//factory.insertPTM(138, new char[]{'K'}, 1);
		//factory.insertPTM(0.9847, new char[]{'N', 'Q'}, 2);
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\20100215_nLCFTMSMS_CLPL_2pmol_1841_MSonly.mzXML";
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\Crosslinked_peps_lib-2_cen.mzXML";
		String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\doug_linkedSinglePeptide\\11_0128_k_dmso_dtt.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		PrecursorMassChecker checker = new PrecursorMassChecker(spectrumLibFile);
		TheoreticalSpectrum.prefixIons = Mass.standardPrefixes;
		TheoreticalSpectrum.suffixIons = Mass.standardSuffixes;
		System.out.println("total library spectrum: " + specList.size());
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.scanNumber > 3600 || s.parentMass < 450){
					continue;					
			}else{
				System.out.println("Query: " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
			}
			List<Peptide> candidates  = factory.getCandidatePeptideByMass(s, 2	, false);
			int candidateCount = 0;
			//s.shiftSpectrumPPM(-100);
			if(candidates.size() > 0){
				for(int j =0; j < candidates.size(); j++){
					Peptide p = candidates.get(j);
					if(checker.matchPrecursorProfile(s.scanNumber, p, 70) < 3){
						//continue;
					}
					System.out.println(s.spectrumName + "\thas peakcount: " + s.getPeaks().size());
					s.windowFilterPeaks(8, 25);
					s.computePeakRank();
					System.out.println("peptide is: " + p.getPeptide());
					TheoreticalSpectrum th = new TheoreticalSpectrum(p);
					double[] stat = th.analyzeAnnotation(s, p.getPeptide(), 0.45);
					System.out.println("Spectrum: " + s.spectrumName+ "\t"  + s.parentMass + " has best match: " +  p + "\t" +  p.getParentmass()
							+ "\t" + stat[0] + "\t" + stat[1] + "\t"
							+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
							+ stat[5]);
					candidateCount++;
				}
			}
			System.out.println(s.spectrumName +"\t" + s.parentMass + "\t" + s.charge + "\thas\t" + candidateCount + " precise single-candidates");
		}		
	}
	
	public static void testGetIsotopePair(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("TR[TN][PG][AY]K[EF][IQ][DS]R");
		List<String> peptides = combPeps.generateAllPeptides();
		addMod(peptides, 42.010564686);
		peptides = combPeps.cutPeptides(peptides);
		//peptides = addMod(peptides);
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(3);
		factory.indexPeptideByParentMass(0.1);
		String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\20100215_nLCFTMSMS_CLPL_2pmol_1841_MSonly.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		System.out.println("total library spectrum: " + specList.size());
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.scanNumber > 3300 || s.parentMass < 450){
				continue;
			}else{
				//reader.isIsotopeCoded(s.scanNumber, 12, 10);
				int charge = reader.getPrecursorCharge(s.scanNumber);
				System.out.println(s.spectrumName + "\t" + s.parentMass + "\t" + charge);
			}	
			
		}		
	}
	
	public static void testMatchTheoreticalMasses(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]LKT[AL]G[FD]V[YA]R");
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]VKL[AL]G[FD]T[YA]R"); 
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("KA[ARDEHLKMFPSTYV]D[ARDEHLKMFPSTYV]ES[ARDEHLKMFPSTYV]LRAK"); 
		List<String> peptides = combPeps.generateAllPeptides();
		//CombinatoryPeptides combPeps2 = new CombinatoryPeptides("QQQTGG");
		CombinatoryPeptides combPeps2 = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		List<String> peptides2 = combPeps2.generateAllPeptides();
		//List<String> peptides2 = Utils.FileIOUtils.createListFromFile("..\\mixture_linked\\lib_albumin_peptides.txt");
		//peptides = appendMod(peptides, 42.010564686);
		//peptides.addAll(peptides2);
		//double linkermass = Math.random()*200;
		//Mass.DSSLINKER_MASS = linkermass;
		System.out.println("we have peptides: " + peptides.size());
		//peptides = addMod(peptides);
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);	
		factory.indexPeptideByParentMass(0.1);
		CandidateSpectrumLibFactory factory2 = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides2);
		factory2.setMinCharge(1);
		factory2.setMaxCharge(1);	
		factory2.indexPeptideByParentMass(0.1);
		//factory.reIndexPeptideByParentMass(0.1);
		//factory.insertPTM(Mass.DSSDANGLE_MASS, new char[]{'K'}, 1);
		//factory.insertPTM(Mass.DSSDANGLE_MASS_D12, new char[]{'K'}, 1);
		//factory.insertPTM(0.9847, new char[]{'N', 'Q'}, 2);
//		factory.crossLinkAllPeptides(3, 3);
		Mass.DSSLINKER_MASS = 138.0680;
		factory.crossLinkAllPeptides(factory2, 2, 4, 4, 4);
//		factory2.crossLinkAllPeptides(3, 3);
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\20100215_nLCFTMSMS_CLPL_2pmol_1841_MSonly.mzXML";
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\Crosslinked_peps_lib-2_cen.mzXML";
		String spectrumLibFile = "..\\mixture_linked\\Sample26masstagCID_101209104646.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		PrecursorMassChecker checker = new PrecursorMassChecker(spectrumLibFile);
		//checker.setSystematicError(-10);
		System.out.println("total library spectrum: " + specList.size());
		List<Peptide> pepList = factory.getAllPeptide();
		long start = (new GregorianCalendar()).getTimeInMillis();
		System.out.println("we have total "+ pepList.size()+ "\tpeptides");
		for(int i = 0; i < pepList.size(); i++){
			Peptide current = pepList.get(i);
			System.out.println("checking : " + current + "\t" + current.getCharge());
			//int matched = checker.matchPeptidePrecursorProfile2(current, 10);
			int[] matched = checker.matchPeptidePrecursorProfilePair(current, 10, (Mass.C13-Mass.C12)*6);
			System.out.println("Peptide " + current +"\t"+current.getCharge()+"\tmatch-stat:\t" + matched[0] + "\t" 
					+ matched[1] + "\t" + matched[2] + "\t" + matched[3] + "\t" + matched[4] + "\t" + matched[5]
					+ "\t" + matched[6]);
		}		
		System.out.println("matching spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	
	public static void testMatchMS3(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]LKT[AL]G[FD]V[YA]R");
		List<String> peptides = combPeps.generateAllPeptides();
		CandidateSpectrumLibFactory factory = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);		
		factory.setMinCharge(1);
		factory.setMaxCharge(1);
		factory.indexPeptideByParentMass(0.5);
		//factory.insertPTM(0.9847, new char[]{'N', 'Q'}, 2);
		factory.crossLinkAllPeptides(4, 6);	
		String training2 = "..\\MSPLib\\Lib\\yeast.msp";
        SpectrumComparator comp2 = SpectrumUtil.getLPeakRankBaseScorer(training2);
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\20100215_nLCFTMSMS_CLPL_2pmol_1841_MSonly.mzXML";
		String spectrumLibFile = "..\\mixture_linked\\Linked_peptide_library\\pepA3_10_0823.mzxml";
		PrecursorMassChecker checker = new PrecursorMassChecker(spectrumLibFile);
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> specList = reader.readAllMS2Spectra();
		reader = null;
		System.out.println("total library spectrum: " + specList.size());
		for(int i = 0; i < specList.size(); i++){
			Spectrum s = specList.get(i);
			if(s.scanNumber > 3600 || s.parentMass < 450 ){
				continue;
			}else{
				System.out.println("Query: " + s.spectrumName + "\t" + s.parentMass + "\t" + s.charge);
			}
			List<Peptide> candidates2  = factory.getCandidatePeptideByMass(s, 2, false);
			System.out.println(s.spectrumName + " " + s.parentMass + " has " + candidates2.size() + " linked-candidates");
			int candidateCount=0;
			if(candidates2.size() > 0){
				//s.shiftSpectrumPPM(-100);
				for(int j =0; j < candidates2.size(); j++){
					s.windowFilterPeaks(8, 25);
					s.computePeakRank();
					LinkedPeptide p = (LinkedPeptide)candidates2.get(j);
					if(checker.matchPrecursorProfile(s.scanNumber, p, 70) < 3){
						//continue;
					}
					TheoreticalSpectrum t = generateLinkedTheoreticalSpectrum(p);
					matchMS2(t	, s, comp2, 0.5);
					candidateCount++;
				}			
			System.out.println(s.spectrumName + "\t" + s.parentMass + "\t" + s.charge + "\thas\t" +  candidateCount + " precise linked-candidates");
			}
		}		
	}
	
	public static TheoreticalSpectrum generateLinkedTheoreticalSpectrum(LinkedPeptide p){
		Peptide p1 = p.peptides[0];
		Peptide p2 = p.peptides[1];
		String[] peps = p.getPeptide().split("--");
		//System.out.println("peptide1 is: " + p1 + "\t" + "peptides2 is: " + p2);
		TheoreticalSpectrum linkedSpect = new TheoreticalSpectrum(p1, p2, p.getCharge(), false); //note linkedSpec parent mass seems not correct in this case
		linkedSpect.p = p;
		//TheoreticalSpectrum linkedSpect = new TheoreticalSpectrum(p, p.getCharge(), false); //note linkedSpec parent mass seems not correct in this case
		//System.out.println("number of peaks genearate for peptide with charge: " + linkedSpect.charge + " is: " + linkedSpect.getPeak().size());
		return linkedSpect;
	}
	
	public static void checkSinglePeptideMass(LinkedPeptide p, Spectrum s, double tolerance){
		//s.filterPeaks(10);
		Peptide p1 = new Peptide(p.peptides[0].getPeptide(),1);
		Peptide p2 = new Peptide(p.peptides[1].getPeptide(),1);
		List<Peak> peakList = s.getPeaks();
		double minDiff = Double.MAX_VALUE, minDiff2 = Double.MAX_VALUE;
		Peak closest = null, closest2 = null;
		for(int i = 0; i < peakList.size(); i++){
			Peak current = peakList.get(i);
			double diff = p1.getParentmass() - current.getMass();
			closest = Math.abs(diff) < Math.abs(minDiff) ? current : closest; 
			minDiff = Math.abs(diff) < Math.abs(minDiff) ? diff : minDiff;
			if(Math.abs(diff) < tolerance){
				System.out.println("Scan " + s.spectrumName + " Precurosr of " 
						+ p1.getPeptide() + "\t" + p1.getParentmass() + " matched to peak " 
						+ current.getIntensity() + "\t"  + current.getMass() + "\t" + current.getRank());
			}
			double diff2 = p2.getParentmass() - current.getMass();
			closest2 = Math.abs(diff2) < Math.abs(minDiff2) ? current : closest2; 
			minDiff2 = Math.abs(diff2) < Math.abs(minDiff2) ? diff2 : minDiff2;
			if(Math.abs(p2.getParentmass() - current.getMass()) < tolerance){
				System.out.println("Scan " + s.spectrumName + " Precurosr of " 
						+ p2.getPeptide() + "\t" + p2.getParentmass() + " matched to peak " 
						+ current.getIntensity() + "\t"  + current.getMass() + "\t" + current.getRank());
			}
		}
		System.out.println("Scan " + s.spectrumName + "\t" + s.parentMass + "\t" + p.getParentmass() +" Precurosr of " 
				+ p1.getPeptide() + "\t" + p1.getParentmass() + " matched to closest-peak " 
				+ closest.getIntensity() + "\t"  + closest.getMass() + "\t" + minDiff);
		System.out.println("Scan " + s.spectrumName + "\t" + s.parentMass + "\t" + p.getParentmass() + " Precurosr of " 
				+ p2.getPeptide() + "\t" + p2.getParentmass() + " matched to closest-peak " 
				+ closest2.getIntensity() + "\t"  + closest2.getMass() + "\t" + minDiff2);
		
	}
	
	public static void matchMS2(TheoreticalSpectrum linkedSpect, Spectrum s, SpectrumComparator comp2, double tolerance){
//		double score = comp2.compare(linkedSpect, s);
//		linkedSpect.analyzeAnnotation(s, linkedSpect.peptide, 0.5);
//		System.out.println(s.spectrumName + "\tscore:\t" + score);
		double score = 0.001;
		if(comp2 != null){
			score = comp2.compare(linkedSpect, s);
		}
		if(score > -100){
			double[] stat = linkedSpect.analyzeMixtureAnnotation(s, "", "", tolerance, false);
			double[] stat2 = MixtureSpectrumComparator.getMatchedIntensity(s, linkedSpect);
			if(!(stat[0] > 0.35 && stat[1] + stat[2] >= 1 && stat[3] + stat[4] >= 1)){
				return;
			}
			linkedSpect.analyzeMixtureAnnotation(s, "", "", tolerance, true);
			System.out.println("Spectrum: " + s.spectrumName+ "\t"  + s.parentMass + "\t" + s.charge+ " has best match: " +  linkedSpect.getPeptide() + "\t" +  linkedSpect.parentMass
				+ "\t" + stat[0] + "\t" + stat[1] + "\t"
				+ stat[2] + "\t" + stat[3] + "\t" + stat[4] + "\t" 
				+ stat[5] + "\t" + stat[6] + "\t" + stat[7] + "\t" + stat[8] + "\t" 
				+ stat[9] + "\t" + stat[10] + "\t"  + stat[13] + "\t" + stat[14] + "\t" 
				+ stat[15] + "\t" + stat[16] + "\t"
				+ stat[11] + "\t" + stat[12] + "\t"
				+ stat2[0] + "\t" + stat2[1] + "\t" + stat2[2] + "\t" + score);
		}
	}
		
	public static List<Peak> getTopMatchedPeak(TheoreticalSpectrum t, Spectrum s, double tolerance){
		List<Peak> pList = SpectrumUtil.getAnnotatedPeak(s, t, tolerance);
		pList = SpectrumUtil.getTopPeaks(2, pList);
		return pList;
	}
	
	//get peaks that add up to parentmass within tolerance
	public static List<Peak> getComplementaryPeaks(Spectrum s, double tolerance){
		List<Peak> peakList = s.getPeaks();
		List<Peak> complements = new ArrayList<Peak>();
		for(int i = 0; i < peakList.size(); i++){
			Peak p1 = peakList.get(i);
			for(int j = i+1; j < peakList.size(); j++){
				Peak p2 = peakList.get(j);
				if(Math.abs(p1.getMass() + p2.getMass() 
						- (s.parentMass*s.charge - (s.charge-2)*Mass.PROTON_MASS)) < tolerance){
					complements.add(p1);
					complements.add(p2);
				}
			}
		}
		System.out.print(s.spectrumName + "\t" + s.parentMass + " complement peaks: ");
		for(int i = 0; i < complements.size(); i++){
			System.out.print(complements.get(i).getMass() + "\t");
		}
		System.out.println();
		return complements;
	}
	
	public static void matchMS3(LinkedPeptide p, Spectrum s, LabelledPeak matchedPeak, SpectrumComparator comp2){
		
	}
	
	public static void analyzeAnnotatedStat(MZXMLReader reader, Map<String, String> annot){
		int count1= 0, count2 = 0, count3 = 0;
		for(Iterator<String> it = annot.keySet().iterator(); it.hasNext();){
			int scan = Integer.parseInt(it.next());
			if(reader.getParser().rap(scan).getHeader().getMsLevel() == 1){
				count1++;
			}
			if(reader.getParser().rap(scan).getHeader().getMsLevel() == 2){
				count2++;
			}
			if(reader.getParser().rap(scan).getHeader().getMsLevel() == 3){
				count3++;
			}			
		}
		System.out.println("Total annotated MS1: " + count1 + " MS2: " + count2 + " MS3: " + count3);
	}
	
	public static void groupMS3Annotation(String spectrumFile, String annotationFile){
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		TreeModel tree = reader.getStructuredSpectrum();		
		System.out.println("Total MS1 gathered: " + tree.getChildCount(tree.getRoot()));
		reader.getSpectrumStat();
		Map<String, String> table = FileIOUtils.createTableFromFile(annotationFile, 1, 2);
		analyzeAnnotatedStat(reader, table);
		List<List<DefaultMutableTreeNode>> MS3Groups = new ArrayList<List<DefaultMutableTreeNode>>();
		List<List<DefaultMutableTreeNode>> MS2Groups = new ArrayList<List<DefaultMutableTreeNode>>();
		List<DefaultMutableTreeNode> MS2s = new ArrayList<DefaultMutableTreeNode>();
		int total = tree.getChildCount(tree.getRoot());
		Object root = tree.getRoot();
		for(int i = 0; i < total; i++){
			DefaultMutableTreeNode ms1 = (DefaultMutableTreeNode)tree.getChild(root, i);
			for(Enumeration<DefaultMutableTreeNode> e = ms1.children(); e.hasMoreElements();){
				MS2s.add(e.nextElement());
			}
		}
		
		List<DefaultMutableTreeNode> toBeRemove, ms3s;
		while(MS2s.size() > 0){
			toBeRemove = new ArrayList<DefaultMutableTreeNode>();
			ms3s = new ArrayList<DefaultMutableTreeNode>();
			DefaultMutableTreeNode node1 = MS2s.get(0);
			Spectrum s1 = (Spectrum)node1.getUserObject();
			toBeRemove.add(node1);
			ms3s.addAll(getChildren(node1));
			for(int j = 1; j < MS2s.size();j++){
				DefaultMutableTreeNode node2 = MS2s.get(j);
				Spectrum s2 = (Spectrum)node2.getUserObject();
				if(Math.abs(s1.parentMass - s2.parentMass) < 0.01
						&& Math.abs(s1.scanNumber - s2.scanNumber) < 100){
					ms3s.addAll(getChildren(node2));
					toBeRemove.add(node2);
				}
			}
			MS2Groups.add(toBeRemove);
			if(ms3s.size() > 0){
				MS3Groups.add(ms3s);
			}
			MS2s.removeAll(toBeRemove);
		}
		System.out.println("There are total of " + MS2Groups.size() + " MS2 groups");
		System.out.println("There are total of " + MS3Groups.size() + " MS3 groups");
		int count = 0;
		for(int i = 0; i < MS3Groups.size(); i++){
			StringBuffer buff = new StringBuffer();
			buff.append("Group " + i + " IDs: ");
			List<DefaultMutableTreeNode> l = MS3Groups.get(i);
			for(int j = 0; j < l.size(); j++){
				Spectrum s = (Spectrum)l.get(j).getUserObject();
				String key = Integer.toString(s.scanNumber);
				//System.out.println("Scan is: " + key);
				if(table.containsKey(key)){
					//System.out.println("found");
					buff.append("Scan: " + key + "\t" + table.get(key) + "\t");
				}
				count++;
			}
			if(buff.length() > 18){
				System.out.println(buff);
			}
		}
		System.out.println("Total processed MS3: " + count);
	}
	
	
	public static void printAnnotationStat(){
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]LKT[AL]G[FD]V[YA]R");
		CombinatoryPeptides combPeps = new CombinatoryPeptides("[YQ]A[EG]V[EL]VKL[AL]G[FD]T[YA]R"); 
		List<String> peptides = combPeps.generateAllPeptides();
		List<String> peptides2 = Utils.FileIOUtils.createListFromFile("..\\mixture_linked\\lib_albumin_peptides.txt");
		peptides = appendMod(peptides, 42.010564686);
		//System.out.println("total peptides: " + peptides.size());
		//CandidateSpectrumLibFactory factory = indexCrossLinkedPeptides(peptides,3,6);
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);	
		factory.indexPeptideByParentMass(0.1);
		CandidateSpectrumLibFactory factory2 = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides2);
		factory2.setMinCharge(1);
		factory2.setMaxCharge(1);	
		factory2.indexPeptideByParentMass(0.1);
		factory.crossLinkAllPeptides(factory2, 3, 3);
//		factory.crossLinkAllPeptides(3, 3);
		factory2.setMinCharge(1);
		factory2.setMaxCharge(3);	
		factory2.indexPeptideByParentMass(0.1);
		//factory2.insertPTM(Mass.DSSDANGLE_MASS, new char[]{'K'}, 1);
		factory2.insertPTM(Mass.DSSLINKER_MASS, new char[]{'K'}, 1);
		String spectrumLibFile = "..\\mixture_linked\\Linked_peptide_library\\lib_ucsd2\\Set9\\Library-clean.mzxml";
		PrecursorMassChecker checker = new PrecursorMassChecker(spectrumLibFile);
		//checker.setSystematicError(-20);
		PreciseCandidatesFactory pFactory = new PreciseCandidatesFactory(factory, checker);
		PreciseCandidatesFactory pFactory2 = new PreciseCandidatesFactory(factory, checker);
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		String annotationFile = "..\\mixture_linked\\MS3_inspect_filtered.txt";
		Map<String, String> table = FileIOUtils.createTableFromFile(annotationFile, 1, 2);
		MSXMLParser parser = reader.getParser();
		for(int i = 1; i < parser.getScanCount(); i++){
			Scan scan = parser.rap(i);
			String key = Integer.toString(i);
			if(scan.getHeader().getMsLevel() == 1){
				continue;
			}
			int index =	reader.getPrevScan(i, scan.getHeader().getMsLevel()-1);
			Spectrum s = reader.getSpectrum(i);
			Spectrum parent = reader.getSpectrum(index);
			List<Peptide> preciseCandidates = null;
			List<Peptide> preciseCandidates2 = null;
			if(scan.getHeader().getMsLevel() == 3){
				preciseCandidates = pFactory.getCandidateByMass(parent, 10);	
				preciseCandidates2 = pFactory2.getCandidateByMass(parent, 10);	
			}else{
				preciseCandidates = pFactory.getCandidateByMass(s, 10);	
				preciseCandidates2 = pFactory2.getCandidateByMass(s, 10);	
			}
			if(table.containsKey(key)){
				String annotation = table.get(key);
				String[] tokens = annotation.split("\\.");
				Peptide p = new Peptide("K",1);
				try{
					p = new Peptide(tokens[1], s.charge);
				}catch(Exception e){
					
				}
				System.out.println("Spectrum: " + "\t" + spectrumLibFile +"\t" + s.scanNumber+ "\t"  
					+ scan.getHeader().getMsLevel() + "\t" + s.parentMass + "\t" + s.charge + "\t" + parent.scanNumber  + "\t"
					+ table.get(key) + "\t" + preciseCandidates.size() + "\t" + preciseCandidates2.size() + "\t" + p.getParentmass() + "\t" 
					+ (p.getParentmass() - s.parentMass)*1000000/p.getParentmass());
			}else{
				System.out.println("Spectrum: " + "\t" + spectrumLibFile +"\t" + s.scanNumber+ "\t"  
					+ scan.getHeader().getMsLevel() + "\t" + s.parentMass + "\t" + s.charge + "\t" + parent.scanNumber  + "\t"
					+ "No-annotation" + "\t" + preciseCandidates.size() + "\t" + preciseCandidates2.size());
			}
			if(preciseCandidates.size() > 0){
				for(int k = 0; k < preciseCandidates.size(); k++){
					Peptide p = preciseCandidates.get(k);
					System.out.println(spectrumLibFile +"\t" + s.scanNumber+ "\t"  
							+ scan.getHeader().getMsLevel() + "\t" + s.parentMass + "\t" + s.charge + "\t" + p + p.getParentmass() + p.getCharge()); 
				}
			}
		}
	}
	
	
	public static List<DefaultMutableTreeNode> getChildren(DefaultMutableTreeNode node){
		List<DefaultMutableTreeNode> children = new ArrayList<DefaultMutableTreeNode>();
		for(Enumeration<DefaultMutableTreeNode> e = node.children(); e.hasMoreElements();){
			children.add(e.nextElement());
		}
		return children;
	}
	
	public static void testGroupMS3Annotation(){
		String spectrumFile = "..\\mixture_linked\\linked_peptide_library\\lib_ucsd2\\Set3\\A1-b.mzxml";
		String annotationFile = "..\\mixture_linked\\MS3_inspect_filtered.txt";		
		groupMS3Annotation(spectrumFile, annotationFile);
	}
	
	public static void generateLibrarySeq(){
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("KA[ARDEHLKMFPSTYV]D[ARDEHLKMFPSTYV]ES[ARDEHLKMFPSTYV]LRAK"); 
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("TALH[ARDEHLKMFPSTYV]K[ARDEHLKMFPSTYV]S[ARDEHLKMFPSTYV]TFR"); 
		CombinatoryPeptides combPeps = new CombinatoryPeptides("A[FILVYW]K[ARDEHLKMFPSTYV][DE]T[ARDEHLKMFPSTYV]FRAK"); 
		List<String> peptides = combPeps.generateAllPeptides();
		for(int i = 0; i < peptides.size(); i++){
			System.out.println(">Peptide: "+peptides.get(i));
			System.out.println(peptides.get(i));
		}
	}
	
	public static void main(String[] args){
		//testGeneratePeptide();
		//testMatchPeptide();
		testMatchPeptideMasses();
		//testMatchSingle();
		//testGetIsotopePair();
		//testMatchTheoreticalMasses();
		//testGroupMS3Annotation();
		//printAnnotationStat();
		//generateLibrarySeq();
	}
}
