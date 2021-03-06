package org.Spectrums;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import org.apache.commons.collections.MultiMap;
import org.apache.commons.collections.map.MultiValueMap;
import Utils.FileIOUtils;


/**
 * Contain method that construct a candidate spectrumlibrary
 * peptides are store with simple indexing by their parent mass
 * @author jian wang
 *
 */
public class CandidateSpectrumLibFactory {
	public Map<Long, Peptide> getPeptideTable() {
		return peptideTable;
	}


	public void setPeptideTable(Map<Long, Peptide> peptideTable) {
		this.peptideTable = peptideTable;
	}

	protected int minCharge = 1;
	protected int maxCharge = 3;
	protected String peptideFile;
	protected String proteinFile;
	protected List<String> peptides;
	protected Map<Long, Peptide> peptideTable;
	protected double massBinWidth;
	protected String[] prefix = {"b"};
	protected String[] suffix = {"y"};
	protected boolean matchCharge = false;
	
	public boolean isMatchCharge() {
		return matchCharge;
	}


	public void setMatchCharge(boolean matchCharge) {
		this.matchCharge = matchCharge;
	}


	public int getMinCharge() {
		return minCharge;
	}


	public void setMinCharge(int minCharge) {
		this.minCharge = minCharge;
	}


	public int getMaxCharge() {
		return maxCharge;
	}


	public void setMaxCharge(int maxCharge) {
		this.maxCharge = maxCharge;
	}


	public String[] getPrefix() {
		return prefix;
	}


	public void setPrefix(String[] prefix) {
		this.prefix = prefix;
	}


	public String[] getSuffix() {
		return suffix;
	}


	public void setSuffix(String[] suffix) {
		this.suffix = suffix;
	}


	protected CandidateSpectrumLibFactory(){
		this.peptides = new ArrayList<String>();
	}
	
	
	public String getPeptideFile() {
		return peptideFile;
	}
	
	public void setPeptideFile(String peptideFile) {
		this.peptideFile = peptideFile;
	}
	
	public String getProteinFile() {
		return proteinFile;
	}
	
	public void setProteinFile(String proteinFile) {
		this.proteinFile = proteinFile;
	}
	
	public List<String> getPeptides() {
		return peptides;
	}
	
	public void setPeptides(List<String> peptides) {
		this.peptides = peptides;
	}
	
	public static CandidateSpectrumLibFactory createFactoryFromPeptide(String peptideFile){
		 CandidateSpectrumLibFactory f = new CandidateSpectrumLibFactory();
		 f.setPeptideFile(peptideFile);
		 f.loadPeptidesFromFile(peptideFile);
		 System.out.println("Done loading peptide");
		 return f;
	}
	
	public static CandidateSpectrumLibFactory createFactoryFromPeptide(List<String> peptides){
		 CandidateSpectrumLibFactory f = new CandidateSpectrumLibFactory();
		 f.setPeptideFile("");
		 f.setPeptides(peptides);
		 return f;
	}
	
	public static CandidateSpectrumLibFactory createFactoryFromProtein(String proteinFile){
		 CandidateSpectrumLibFactory f = new CandidateSpectrumLibFactory();
		 f.setPeptideFile(proteinFile);
		 //need to use Digester to digest proteins before use
		 return f;
	}
	
	
	public void loadPeptidesFromFile(String file){
		this.setPeptideFile(file);
//		Map<String, String> table = FileIOUtils.createTableFromFile(file, 0, 0);
		List<String> peps = FileIOUtils.createListFromFile(file);
		this.peptides.addAll(peps);
	}
	
	
	public void indexPeptideByParentMass(){
		indexPeptideByParentMass(0.5);
	}
	
	public void indexPeptideByParentMass(double massBinWidth){
		Peptide p;
		Long massIndex;
		Collection<Peptide> pepList;
		this.peptideTable = new MultiValueMap();
		this.massBinWidth = massBinWidth;
		//System.out.println("we start with " + this.peptides.size() + " peptides");
		this.peptideTable.clear();
		for(int i = 0, size = this.peptides.size(); i < size; i++){
			for(int charge = this.minCharge; 
				charge <= this.maxCharge;
				charge++){
				String str = this.peptides.get(i);
				if(str.startsWith("r")){
					p = new Peptide(str.substring(1), charge);
//					p.setDecoy(true);
				}else{
					//p  = new Peptide(this.peptides.get(i), charge);
					p  = new Peptide(this.peptides.get(i), charge);
				}
				//System.out.println("peptide is: " + p + "\t" + p.getParentmass() + "\t" + p.getCharge());
				long index = Math.round((float)p.getParentmass()/massBinWidth);
				massIndex = new Long(index);
				//System.out.println("key is " + mass);
				this.peptideTable.put(index, p);
			}
		}
		System.out.println("indexed " + this.peptideTable.values().size() + " peptides");
	}
	
	public void indexPeptideByParentMass(double massBinWidth, List<Peptide> pepList){
		Peptide p;
		Long massIndex;
		this.peptideTable = new MultiValueMap();
		this.massBinWidth = massBinWidth;
		//System.out.println("we start with " + this.peptides.size() + " peptides");
		this.peptideTable.clear();
		for(int i = 0, size = pepList.size(); i < size; i++){
				p = pepList.get(i);
				long index = Math.round((float)p.getParentmass()/massBinWidth);
				massIndex = new Long(index);
				//System.out.println("key is " + mass);
				this.peptideTable.put(index, p);

		}
		System.out.println("indexed " + this.peptideTable.values().size() + " peptides");
	}
	
	
	//be-careful of this method
	public void reIndexPeptideByParentMass(double massBinWidth){
		Long massIndex;
		List<Peptide> pepList;
		HashMap<Long, List<Peptide>> newTable = new HashMap<Long, List<Peptide>>();
		this.massBinWidth = massBinWidth;	
		//System.out.println("we start with " + this.peptides.size() + " peptides");
		List<Peptide> allpeptides = new ArrayList<Peptide>();
		allpeptides.addAll(((MultiMap)this.peptideTable).values());
		this.peptideTable.clear();
		for(int i = 0; i < allpeptides.size(); i++){
				Peptide p = allpeptides.get(i);
				//p.insertPTM(p.getPeptide().indexOf('K')+1, Mass.DSSDANGLE_MASS);
				long index = Math.round((float)p.getParentmass()/massBinWidth);
				massIndex = new Long(index);
				//System.out.println("key is " + mass);
				this.peptideTable.put(massIndex, p);
		}
		System.out.println("indexed " + this.peptideTable.values().size() + " peptides");
	}
	
	public void insertPTM(double mass, char[] residues, int maxPTM){
		List<Peptide> allpeptides = new ArrayList<Peptide>();
		allpeptides.addAll(((MultiMap)this.peptideTable).values());
		System.out.println("we start with " + allpeptides.size());
		List<Peptide> toBeAdded = Peptide.insertPTM(allpeptides, mass, residues, maxPTM);
		allpeptides.addAll(toBeAdded);
		System.out.println("we end up with " + allpeptides.size());
		reIndexPeptideByParentMass(this.massBinWidth, allpeptides);
	}
	
	public void insertPTM(double mass, int position, int maxPTM){
		List<Peptide> allpeptides = new ArrayList<Peptide>();
		allpeptides.addAll(((MultiMap)this.peptideTable).values());
		System.out.println("we start with " + allpeptides.size());
		List<Peptide> toBeAdded = Peptide.insertPTM(allpeptides, mass, position, maxPTM);
		allpeptides.addAll(toBeAdded);
		System.out.println("we end up with " + allpeptides.size());
		reIndexPeptideByParentMass(this.massBinWidth, allpeptides);
	}
	
	public void crossLinkAllPeptides(int minCharge, int maxCharge){
//		List<Peptide> allpeptides = new ArrayList<Peptide>();
//		allpeptides.addAll(((MultiMap)this.peptideTable).values());
//		System.out.println("we start with " + allpeptides.size() + " peptides");
//		List<Peptide> toBeAdded = new ArrayList<Peptide>();
//		for(int i = 0; i < allpeptides.size(); i++){
//			Peptide p1 = allpeptides.get(i);
//			List<Integer> positions = p1.getLysPositions();
//			for(int j = i+1; j < allpeptides.size(); j++){
//				Peptide p2 = allpeptides.get(j);
//				List<Integer> positions2 = p2.getLysPositions();
//				for(int c = minCharge; c <= maxCharge; c++){
//					for(int m = 0; m < positions.size(); m++){
//						for(int n = 0; n < positions2.size(); n++){
//							LinkedPeptide linked = new LinkedPeptide(p1, p2, c, positions.get(m)+1, positions2.get(n)+1);
//							toBeAdded.add(linked);
//						}
//					}
//				}
//			}
//		}
		List<Peptide> allpeptides = this.getAllPeptide();
		this.minCharge = minCharge;
		this.maxCharge = maxCharge;
		crossLinkAllPeptides(allpeptides, allpeptides);
//		reIndexPeptideByParentMass(this.massBinWidth, toBeAdded);
		
	}
	
	public void crossLinkAllPeptides(CandidateSpectrumLibFactory lib2, int minCharge, int maxCharge){
		List<Peptide> allpeptides = this.getAllPeptide();
		List<Peptide> allpeptides2 = lib2.getAllPeptide();
		this.minCharge = minCharge;
		this.maxCharge = maxCharge;
		crossLinkAllPeptides(allpeptides, allpeptides2);
	}
	
	public void crossLinkAllPeptides(CandidateSpectrumLibFactory lib2, int minCharge, int maxCharge, int position1, int position2){
		List<Peptide> allpeptides = this.getAllPeptide();
		List<Peptide> allpeptides2 = lib2.getAllPeptide();
		this.minCharge = minCharge;
		this.maxCharge = maxCharge;
		crossLinkAllPeptides(allpeptides, allpeptides2, position1, position2);
	}
	
	public void crossLinkAllPeptides(CandidateSpectrumLibFactory lib2, int minCharge, int maxCharge, CrossLinker linker){
		List<Peptide> allpeptides = this.getAllPeptide();
		List<Peptide> allpeptides2 = lib2.getAllPeptide();
		this.minCharge = minCharge;
		this.maxCharge = maxCharge;
		crossLinkAllPeptides(allpeptides, allpeptides2, linker);
	}
	
	public List<Peptide> crossLinkAllPeptides(List<Peptide> list1, List<Peptide>list2){
		List<Peptide> toBeAdded = new ArrayList<Peptide>();
		for(int i = 0; i < list1.size(); i++){
			Peptide p1 = list1.get(i);
			List<Integer> positions = p1.getLysPositions();
			for(int j = 0; j < list2.size(); j++){
				Peptide p2 = list2.get(j);
				if(p1.getPeptide().equals(p2.getPeptide())){
					continue;
				}
				List<Integer> positions2 = p2.getLysPositions();
				for(int c = minCharge; c <= maxCharge; c++){
					for(int m = 0; m < positions.size(); m++){
						for(int n = 0; n < positions2.size(); n++){
							LinkedPeptide linked = null;
							if(p2.getPeptide().equals("Z")){
								linked = new LinkedPeptide(p1, p2, c, 0, 1);
							}else{
								linked = new LinkedPeptide(p1, p2, c, positions.get(m)+1, positions2.get(n)+1);
							}
							toBeAdded.add(linked);
						}
					}
				}
			}
		}
		reIndexPeptideByParentMass(this.massBinWidth, toBeAdded);
		return toBeAdded;
	}
	
	public List<Peptide> crossLinkAllPeptides(List<Peptide> list1, List<Peptide>list2, CrossLinker linker){
		List<Peptide> toBeAdded = new ArrayList<Peptide>();
		for(int i = 0; i < list1.size(); i++){
			Peptide p1 = list1.get(i);
			List<Integer> positions = linker.getLinkerPositions1(p1.getPeptide());
			for(int j = 0; j < list2.size(); j++){
				Peptide p2 = list2.get(j);
				if(p1 == p2){
					continue;
				}
				List<Integer> positions2 = linker.getLinkerPositions2(p2.getPeptide());
				for(int c = minCharge; c <= maxCharge; c++){
					for(int m = 0; m < positions.size(); m++){
						for(int n = 0; n < positions2.size(); n++){
							LinkedPeptide linked = new LinkedPeptide(p1, p2, c, positions.get(m)+1, positions2.get(n)+1);
							System.out.println(linked +"\t" + linked.getParentmass());
							//System.out.println("number of linked peptides: " + toBeAdded.size());
							toBeAdded.add(linked);
						}
					}
				}
			}
		}
		reIndexPeptideByParentMass(this.massBinWidth, toBeAdded);
		return toBeAdded;
	}
	
	public static List<Peptide> getCrossLinkPeptides(List<Peptide> list1, List<Peptide>list2, CrossLinker linker, int minCharge, int maxCharge){
		List<Peptide> toBeAdded = new ArrayList<Peptide>();
		for(int j = 0; j < list2.size(); j++){
			toBeAdded.addAll(getCrossLinkPeptides(list1, list2.get(j), linker, minCharge, maxCharge));
		}
		return toBeAdded;
	}
	
	public static List<Peptide> getCrossLinkPeptides(List<Peptide> list1, Peptide p2, CrossLinker linker, int minCharge, int maxCharge){
		List<Peptide> toBeAdded = new ArrayList<Peptide>();
		for(int i = 0; i < list1.size(); i++){
			Peptide p1 = list1.get(i);
			List<Integer> positions = linker.getLinkerPositions1(p1.getPeptide());
				if(p1 == p2){
					continue;
				}
				List<Integer> positions2 = linker.getLinkerPositions2(p2.getPeptide());
				for(int c = minCharge; c <= maxCharge; c++){
					for(int m = 0; m < positions.size(); m++){
						for(int n = 0; n < positions2.size(); n++){
							LinkedPeptide linked = new LinkedPeptide(p1, p2, c, positions.get(m)+1, positions2.get(n)+1);
							//System.out.println(p1);
							//System.out.println("number of linked peptides: " + toBeAdded.size());
							toBeAdded.add(linked);
						}
					}
				}
		}
		return toBeAdded;
	}
	
	public static List<Peptide> getCrossLinkPeptides(Peptide p1, Peptide p2, CrossLinker linker, int minCharge, int maxCharge){
		List<Peptide> toBeAdded = new ArrayList<Peptide>();
		List<Integer> positions = linker.getLinkerPositions1(p1.getPeptide());
		List<Integer> positions2 = linker.getLinkerPositions2(p2.getPeptide());
		for(int c = minCharge; c <= maxCharge; c++){
			for(int m = 0; m < positions.size(); m++){
				for(int n = 0; n < positions2.size(); n++){
					LinkedPeptide linked = new LinkedPeptide(p1, p2, c, positions.get(m)+1, positions2.get(n)+1, linker.getLinkerMassOffSet());
					//System.out.println(p1);
					//System.out.println("number of linked peptides: " + toBeAdded.size());
					toBeAdded.add(linked);
				}
			}
		}
		return toBeAdded;
	}
	
	public void internalLinkPeptides(int minCharge, int maxCharge, CrossLinker linker){
		List<Peptide> allpeptides = this.getAllPeptide();
		this.minCharge = minCharge;
		this.maxCharge = maxCharge;
		internalLinkPeptides(allpeptides, linker);
	}
	
	public List<Peptide> internalLinkPeptides(List<Peptide> list1, CrossLinker linker){
		List<Peptide> toBeAdded = new ArrayList<Peptide>();
		for(int i = 0; i < list1.size(); i++){
			Peptide p1 = list1.get(i);
			List<Integer> positions = linker.getLinkerPositions1(p1.getPeptide());
			List<Integer> positions2 = linker.getLinkerPositions2(p1.getPeptide());
				for(int c = minCharge; c <= maxCharge; c++){
					for(int m = 0; m < positions.size(); m++){
						for(int n = 0; n < positions2.size(); n++){
							if(positions.get(m) < positions2.get(n)){
								InternalLinkedPeptide iLinked = new InternalLinkedPeptide(p1, c, positions.get(m)+1, positions2.get(n)+1);
								System.out.println("internal linked peptide: " + iLinked + "\t" + iLinked.getCharge() +  "\t" + iLinked.getParentmass());
								toBeAdded.add(iLinked);
							}
						}
					}
				}
			}
		reIndexPeptideByParentMass(this.massBinWidth, toBeAdded);
		return toBeAdded;
	}
	
	
	public List<Peptide> crossLinkAllPeptides(List<Peptide> list1, List<Peptide>list2, int position1, int position2){
		List<Peptide> toBeAdded = new ArrayList<Peptide>();
		for(int i = 0; i < list1.size(); i++){
			Peptide p1 = list1.get(i);
			List<Integer> positions = p1.getLysPositions();
			for(int j = 0; j < list2.size(); j++){
				Peptide p2 = list2.get(j);
				if(p1.getPeptide().equals(p2.getPeptide())){
					continue;
				}
				List<Integer> positions2 = p2.getLysPositions();
				for(int c = minCharge; c <= maxCharge; c++){
					LinkedPeptide linked = new LinkedPeptide(p1, p2, c, position1, position2);
					toBeAdded.add(linked);
				}
			}
		}
		reIndexPeptideByParentMass(this.massBinWidth, toBeAdded);
		return toBeAdded;
	}
	

	
	public void reIndexPeptideByParentMass(double massBinWidth, List<Peptide> allpeptides){
		Long massIndex;
		List<Peptide> pepList;
		this.massBinWidth = massBinWidth;
		this.peptideTable.clear();
		for(int i = 0; i < allpeptides.size(); i++){
				Peptide p = allpeptides.get(i);
				long index = Math.round((float)p.getParentmass()/massBinWidth);
				massIndex = new Long(index);
				//System.out.println("peptide  is " + p.toString() + " mass: " + p.getParentmass() + "\t" + p.getCharge());
				this.peptideTable.put(massIndex, p);
		}
		System.out.println("indexed " + this.peptideTable.values().size() + " peptides");
	}
	
	public void indexLinkedPeptideByParentMass(double massBinWidth){
		Peptide p;
		Long massIndex;
		List<Peptide> pepList;
		this.peptideTable = new MultiValueMap();
		this.massBinWidth = massBinWidth;
		//System.out.println("we start with " + this.peptides.size() + " peptides");
		for(int i = 0, size = this.peptides.size(); i < size; i++){
			for(int charge = this.minCharge; 
				charge <= this.maxCharge;
				charge++){
				String str = this.peptides.get(i);
				if(str.startsWith("r")){
					p = new Peptide(str.substring(1), charge);
					p.setDecoy(true);
				}else{
					p  = new LinkedPeptide(this.peptides.get(i), charge);
				}
				System.out.println("peptide: " + p.getPeptide() + "\t" + p.getParentmass() + "\t" + p.getCharge());
				long index = Math.round((float)p.getParentmass()/massBinWidth);
				massIndex = new Long(index);
				this.peptideTable.put(massIndex, p);
			}
		}
		System.out.println("indexed " + this.peptideTable.values().size() + " peptides");
	}

	public List<Peptide> getCandidateByMass(double parentMass, double tolerance){
		//this.matchCharge = false;
		return getCandidateByMass(parentMass, 0, tolerance);
	}
	public List<Peptide> getCandidateByMass(double parentMass, int charge, double tolerance){
		//System.out.println("mass is: " + parentMass);
		long pm = Math.round(parentMass/this.massBinWidth);
		long tol = Math.round(tolerance/this.massBinWidth);
		//System.out.println("massindex is: " + pm);
		//System.out.println("tolerance is: " + tol);
		long indexWidth = tol < 0 ? 1 : 2*tol;
		Collection<Peptide> subList; 
		List<Peptide> pepList = new ArrayList<Peptide>();
		Peptide p;
		//System.out.println("left index: " + (pm-indexWidth) + " right index: " + (pm + indexWidth));
		for(Long left = new Long(pm-indexWidth), right = new Long(pm + indexWidth); 
			left <= right; left++ ){
			//System.out.println("checking table index: " + left);
			if(this.peptideTable.containsKey(left)){
				subList = ((MultiValueMap)this.peptideTable).getCollection(left);
				//System.out.println("entry has size: " + subList.size());
				for(Iterator<Peptide> it = subList.iterator(); it.hasNext();){
					p = it.next();
					if(Math.abs(p.getParentmass() - parentMass) < tolerance){
						if(!this.matchCharge || charge == p.getCharge()){
							pepList.add(p);
						}
					}
				}
			}
		}
		return pepList;
	}
	
//	public List<Peptide> getCandidateByMass(double parentMass, double tolerance){
//		//System.out.println("mass is: " + parentMass);
//		long pm = Math.round(parentMass/this.massBinWidth);
//		long tol = Math.round(tolerance/this.massBinWidth);
//		//System.out.println("massindex is: " + pm);
//		//System.out.println("tolerance is: " + tol);
//		long indexWidth = tol < 0 ? 1 : 2*tol;
//		Collection<Peptide> subList; 
//		List<Peptide> pepList = new ArrayList<Peptide>();
//		Peptide p;
//		//System.out.println("left index: " + (pm-indexWidth) + " right index: " + (pm + indexWidth));
//		for(Long left = new Long(pm-indexWidth), right = new Long(pm + indexWidth); 
//			left <= right; left++ ){
//			//System.out.println("checking table index: " + left);
//			if(this.peptideTable.containsKey(left)){
//				subList = ((MultiValueMap)this.peptideTable).getCollection(left);
//				//System.out.println("entry has size: " + subList.size());
//				for(Iterator<Peptide> it = subList.iterator(); it.hasNext();){
//					p = it.next();
//					if(Math.abs(p.getParentmass() - parentMass) < tolerance){
//						pepList.add(p);
//					}
//				}
//			}
//		}
//		return pepList;
//	}
	
	public List<Peptide> getAllPeptide(){
		List<Peptide> pepList = new ArrayList<Peptide>();
		pepList.addAll(((MultiMap)this.peptideTable).values());
		return pepList;
	}
	
	public Peptide getClosestCandidate(double parentMass){
		double min = 1000000000;
		Peptide best = null;
		for(Iterator<Peptide> it = this.getAllPeptide().iterator(); it.hasNext();){
			Peptide current = it.next();
			//double diff = Math.abs(current.getParentmass() - parentMass);
			double diff = Math.abs((current.getParentmass() - parentMass))*1000000/parentMass;
			best = min < diff ? best : current;
			min = min < diff ? min : diff;
		}
		System.out.println("closest matched peptide: " + best.getPeptide() + " with error: " +  min);
		return best;
	}
	
	public SpectrumLib createCandidateSpectrumLib(Spectrum s, double pmTolerance){
		return createCandidateSpectrumLib(s, pmTolerance, true);
	}
	
	public SpectrumLib createCandidateSpectrumLib(Spectrum s, double pmTolerance, boolean isChargeCertain){
		if(this.peptides.size() < 1){
			this.loadPeptidesFromFile(this.peptideFile);
		}
		Map<String,List<Spectrum>> table = new HashMap<String,List<Spectrum>>();
		Iterator<String> it = peptides.iterator();
		Peptide p;
		TheoreticalSpectrum t;
		List<Spectrum> l;
		int startCharge, limitCharge;
		//if we trust the charge of MS2, use the charge information to
		//reduce the candidate list, otherwise we iterate through all 
		//possible charges
		if(isChargeCertain){
			startCharge = s.charge;
			limitCharge = s.charge;
		}else{
			startCharge = this.minCharge;
			limitCharge = this.maxCharge;
		}
	
		String current;
		while(it.hasNext()){
			current = it.next();
			for(int c = startCharge; c <= limitCharge; c++){
				p = new Peptide(current + "." + c);
				if(Math.abs(p.getParentmass() - s.parentMass) < pmTolerance){
					t = new TheoreticalSpectrum(p);
					l = new ArrayList<Spectrum>();
					l.add(t);
					table.put(t.peptide, l);
				}
			}
		}
		return new SpectrumLib(table);
	}
	
	public SpectrumLib createCandidateSpectrumLibX(Spectrum s, double pmTolerance, boolean isDecoy){
		if(this.peptides.size() < 1){
			this.loadPeptidesFromFile(this.peptideFile);
		}
		if(this.peptideTable == null){
			this.indexPeptideByParentMass();
		}
		List<Peptide> candidates = this.getCandidateByMass(s.parentMass, s.charge, pmTolerance);
		//List<Peptide> candidates2 = this.getCandidateByMass(s.parentMass-57, s.charge, pmTolerance);
		//List<Peptide> candMod1 = Peptide.insertPTM(candidates2, -57, 'C', 1);
		//candidates.addAll(candMod1);
		System.out.println("Scan: " + s.scanNumber + "\t"+ s.charge + "\tnumber of candidates: " + candidates.size());
		if(isDecoy){
			return createPlusDecoyFromPeptides(candidates, this.prefix, this.suffix);
		}else{
			return createLibFromPeptides(candidates, this.prefix, this.suffix);
			//return createLazyLibFromPeptides(candidates, this.prefix, this.suffix);
		}
	}


	public List<Peptide> getCandidatePeptideByMass(Spectrum s, double pmTolerance, boolean isDecoy){
		if(this.peptides.size() < 1){
			this.loadPeptidesFromFile(this.peptideFile);
		}
		if(this.peptideTable == null){
			this.indexPeptideByParentMass();
		}
		List<Peptide> candidates = this.getCandidateByMass(s.parentMass, s.charge, pmTolerance);
		return candidates;
	}
	
	
	public SpectrumLib createLibFromPeptides(List<Peptide> pList, String[] prefixIons, String[] suffixIons){
		Map<String,List<Spectrum>> table = new HashMap<String,List<Spectrum>>();
		Iterator<Peptide> it = pList.iterator();
		Peptide current;
		List<Spectrum> l;
		TheoreticalSpectrum t;
		LabelledPeakFactory.setPoolMode(false);
		LabelledPeakFactory.resetFactory();
		while(it.hasNext()){
			current = it.next();
			//System.out.println("candidate is: " + current);
			t = new TheoreticalSpectrum(current, prefixIons, suffixIons);
			l = new ArrayList<Spectrum>();
			l.add(t);
			table.put(t.p.toString(), l);
		}
		return new SpectrumLib(table);
	}
	
	public SpectrumLib createLazyLibFromPeptides(List<Peptide> pList, String[] prefixIons, String[] suffixIons){
		System.out.println("creating lazy lib");
		Map<String,List<Spectrum>> table = new HashMap<String,List<Spectrum>>();
		Iterator<Peptide> it = pList.iterator();
		Peptide current;
		List<Spectrum> l;
		TheoreticalSpectrum t;
		while(it.hasNext()){
			current = it.next();
			t = new LazyEvaluatedSpectrum(current);
			l = new ArrayList<Spectrum>();
			l.add(t);
			table.put(t.peptide, l);
		}
		return new SpectrumLib(table);
	}

	public SpectrumLib createDecoyFromPeptides(List<Peptide> pList, String[] prefixIons, String[] suffixIons){
		Map<String,List<Spectrum>> table = new HashMap<String,List<Spectrum>>();
		Iterator<Peptide> it = pList.iterator();
		Peptide current;
		List<Spectrum> l;
		TheoreticalSpectrum t;
		while(it.hasNext()){
			current = it.next();
			current = current.reverse();
			t = new TheoreticalSpectrum(current, prefixIons, suffixIons);
			l = new ArrayList<Spectrum>();
			l.add(t);
			table.put(t.peptide, l);
		}
		return new SpectrumLib(table);
	}
	
	public SpectrumLib createPlusDecoyFromPeptides(List<Peptide> pList, String[] prefixIons, String[] suffixIons){
		Map<String,List<Spectrum>> table = new HashMap<String,List<Spectrum>>();
		Iterator<Peptide> it = pList.iterator();
		Peptide current;
		Peptide reverse;
		List<Spectrum> l, lr;
		TheoreticalSpectrum t, rt;
		while(it.hasNext()){
			current = it.next();
			reverse = current.reverse();
			reverse.setDecoy(true);
			t = new TheoreticalSpectrum(current, prefixIons, suffixIons);
			rt = new TheoreticalSpectrum(reverse, prefixIons, suffixIons);
			l = new ArrayList<Spectrum>();
			lr = new ArrayList<Spectrum>();
			l.add(t);
			lr.add(rt);
			table.put(rt.peptide, lr);
			table.put(t.peptide, l);
		}
		return new SpectrumLib(table);
	}
	
	public static void testCreatedCandidates(){
		String file = "..\\mixture_linked\\Ecoli_allpeptides.txt";
		CandidateSpectrumLibFactory f = 
			CandidateSpectrumLibFactory.createFactoryFromPeptide(file);
		f.loadPeptidesFromFile(file);
		f.indexPeptideByParentMass();
		Spectrum s = new Spectrum();
		s.charge = 2;
		s.parentMass = 1000.0;
		s.peptide = "KVIITAPAK.2";
		long start = (new GregorianCalendar()).getTimeInMillis();
		SpectrumLib lib = f.createCandidateSpectrumLib(s, 5, false);
		System.out.println("Got " + lib.getSpectrumList().size() + " candidates total");
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
		start = (new GregorianCalendar()).getTimeInMillis();
		SpectrumLib lib2 = f.createCandidateSpectrumLibX(s, 5, false);
		System.out.println("Got " + lib2.getSpectrumList().size() + " candidates total");
		System.out.println("Running for: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void main(String[] args){
		testCreatedCandidates();
	}
	
}
