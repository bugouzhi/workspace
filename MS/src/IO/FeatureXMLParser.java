/**
 * Parse and analyze output from feature finder of openMS
 */
package IO;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Iterator;

import javax.xml.parsers.*;

import org.Spectrums.CandidateSpectrumLibFactory;
import org.Spectrums.CombinatoryPeptides;
import org.Spectrums.CrossLinker;
import org.Spectrums.LinkedPeptide;
import org.Spectrums.MSFeature;
import org.Spectrums.MSFeatureAcq;
import org.Spectrums.Mass;
import org.Spectrums.Peptide;
import org.Spectrums.PrecursorMassChecker;
import org.Spectrums.Spectrum;
import org.Spectrums.SpectrumUtil;
import org.apache.commons.collections.map.MultiValueMap;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import Utils.FileIOUtils;


public class FeatureXMLParser {
	private Document dom;
	private boolean mapMS2Yet=false;
	private List<MSFeature> featureList;
	boolean verbose = true;
	
	public Document getDom() {
		return dom;
	}

	public void setDom(Document dom) {
		this.dom = dom;
	}

	public List<MSFeature> getFeatureList() {
		return featureList;
	}

	public void setFeatureList(List<MSFeature> featureList) {
		this.featureList = featureList;
	}
	
	//create a parser for xml file
	//if input string is folder, parse all files and extrac all features
	public FeatureXMLParser(String filename){
		File path = new File(filename);
		this.featureList = new ArrayList<MSFeature>();
		if(path.isDirectory()){
			File[] files = path.listFiles();
			for(int i = 0; i < files.length; i++){
				String filepath = files[i].getAbsolutePath();
				FeatureXMLParser parser = new FeatureXMLParser(filepath);
				this.featureList.addAll(parser.getFeatureList());
			}
		}else if(path.exists()){
			parseXmlFile(filename);
			parseDocument();
		}
	}
	
	public FeatureXMLParser(String filename, String filter){
		File path = new File(filename);
		this.featureList = new ArrayList<MSFeature>();
		if(path.isDirectory()){
			File[] files = path.listFiles();
			for(int i = 0; i < files.length; i++){
				if(files[i].getName().matches(filter)){
					String filepath = files[i].getAbsolutePath();
					FeatureXMLParser parser = new FeatureXMLParser(filepath);
					this.featureList.addAll(parser.getFeatureList());
				}
			}
		}else if(path.exists()){
			parseXmlFile(filename);
			parseDocument();
		}
	}
	
	private void parseXmlFile(String filename){
		//get the factory
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		
		try {

			//Using factory get an instance of document builder
			DocumentBuilder db = dbf.newDocumentBuilder();

			//parse using builder to get DOM representation of the XML file
			dom = db.parse(filename);
			//dom = db.parse("..//mixture_linked//employee.xml");


		}catch(ParserConfigurationException pce) {
			pce.printStackTrace();
		}catch(SAXException se) {
			se.printStackTrace();
		}catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	private void parseDocument(){
		//get the root element
		Element docEle = dom.getDocumentElement();
		//System.out.println("encoding is: " + dom.getInputEncoding());
		//get a nodelist of 
		
		NodeList nl = docEle.getElementsByTagName("feature");
		System.out.println("There are total of features: " + nl.getLength());
		this.featureList = new ArrayList();
		if(nl != null && nl.getLength() > 0) {
			int count = 0;
			for(int i = 0 ; i < nl.getLength(); i++) {
				//get the employee element
				Element el = (Element)nl.item(i);
				//System.out.println("feature id is: " + el.getAttribute("id"));			
				Node parent = el.getParentNode();
				if(parent.getNodeName().equals("featureList")){
					MSFeature feature = new MSFeatureAcq();
					double[] timeSpan = parseHullPoint(el);
					double[] location = getPositionValue(el, "position");
					feature.setId(el.getAttribute("id"));
					feature.setIntensity(getDoubleValue(el, "intensity"));
					feature.setMz(location[0]);
					feature.setQuality(getDoubleValue(el, "overallquality"));
					feature.setCharge(getIntValue(el, "charge"));
					//feature.setFeatureList(parseFeature(el));
					feature.setMinRT(timeSpan[0]);
					feature.setMaxRT(timeSpan[1]);
					feature.setRt(0.5*(feature.getMaxRT()-feature.getMinRT())+feature.getMinRT());
					feature.setMz(location[1]);
//					System.out.print(el.getAttribute("id") + "\t" 
//							+ getTextValue(el, "intensity") + "\t"
//							+ getTextValue(el, "overallquality") + "\t"
//							+ location[0] + "\t" + location[1] + "\t" 
//							+ getTextValue(el, "charge") + "\t"
//							+ timeSpan[0] +"\t" + timeSpan[1] +"\t"
//							+ "\n");
					featureList.add(feature);
					count++;
				}
				
			}
			System.out.println("there are total of " + count + " main features");
//			for(int i =0; i < featureList.size(); i++){
//				System.out.println(featureList.get(i));
//				
//			}
		}
	}
	
	private MultiValueMap parseFeature(Element el){
		MultiValueMap map = new MultiValueMap();
		List<Element> hullPoints = getTopLevelHullPoints(el);
		System.out.println("extracted hullpoints: " + hullPoints.size());
		for(int i = 0; i < hullPoints.size(); i++){
				Element hullpoint = hullPoints.get(i);
				double[] positions = parseHullPoint(hullpoint);
				map.put(positions[0], positions[1]);
		}
		return map;
	}
	/**
	 * Only extract top level hull points, ignore those in subordinate 
	 * @param el
	 * @return
	 */
	private List<Element> getTopLevelHullPoints(Element el){
		NodeList nl = el.getElementsByTagName("hullpoint");
		List<Element> hullPoints = new ArrayList();
		for(int i = 0; i < nl.getLength(); i++){
			Element hullpoint = (Element)nl.item(i);
			if(hullpoint.getParentNode().getParentNode()
					.getParentNode().getNodeName().equals("featureList")){
				hullPoints.add(hullpoint);
			}
		}
		return hullPoints;
	}
	
	//parse the hullpoint object, extract begin and end elution time of this peptide feature
	private double[] parseHullPoint(Element hullpoint){
		double[] positions = new double[2];
		double beginTime = 0;
		double endTime = 0;
		NodeList nl = hullpoint.getElementsByTagName("convexhull");
		//System.out.println("convexhull size: " + nl.getLength());
		for(int i = 0; i < nl.getLength(); i++){
			Element hullLine = (Element)nl.item(i);
			NodeList pts = hullLine.getElementsByTagName("pt");
			//System.out.println("number of points in hull " + hullLine.getElementsByTagName("pt").getLength());
			for(int j = 0; j < pts.getLength(); j++){
				//Element pt = (Element)pts.item(j);
				//System.out.println("point " + pt.getAttribute("x") + "\t" + pt.getAttribute("y"));
			}
			if(i==0){ //want to get widest time-span, presumably it is the first window width
				beginTime = Double.parseDouble(((Element)pts.item(0)).getAttribute("x"));
				endTime = Double.parseDouble(((Element)pts.item(1)).getAttribute("x"));
			}
		}
		positions[0] = beginTime;
		positions[1] = endTime;
		return positions;
	}
	
	private double[] parsePosition(Element feature){
		double[] positions = new double[2];
		NodeList nl = feature.getElementsByTagName("position");
		positions[0] = Double.parseDouble(nl.item(0).getFirstChild().getNodeValue());
		positions[1] = Double.parseDouble(nl.item(1).getFirstChild().getNodeValue());
		System.out.println("position is: " + positions[0] + "\t" + positions[1]);
		return positions;
	}
	
	private String getTextValue(Element ele, String tagName) {
		String textVal = null;
		NodeList nl = ele.getElementsByTagName(tagName);
		if(nl != null && nl.getLength() > 0) {
			Element el = (Element)nl.item(0);
			textVal = el.getFirstChild().getNodeValue();
		}

		return textVal;
	}
	
	private double[] getPositionValue(Element ele, String tagName) {
		String textVal = null, textVal2 = null;
		NodeList nl = ele.getElementsByTagName(tagName);
		if(nl != null && nl.getLength() > 0) {
			Element el = (Element)nl.item(0);
			textVal = el.getFirstChild().getNodeValue();
			Element el2 = (Element)nl.item(1);
			textVal2 = el2.getFirstChild().getNodeValue();
		}

		return new double[]{Double.parseDouble(textVal), 
				Double.parseDouble(textVal2)};
	}
	
	private int getIntValue(Element ele, String tagName) {
		return Integer.parseInt(getTextValue(ele,tagName));
	}

	private double getDoubleValue(Element ele, String tagName) {
		//in production application you would catch the exception
		return Double.parseDouble(getTextValue(ele,tagName));
	}
	
	public List<MSFeature[]>getFeaturePair(double expectedMassDiff, double massTolerance, double maxRTDiff, int mode){
		List<MSFeature[]> pairedFeatures = new ArrayList<MSFeature[]>();
		Map<String, MSFeature> paired = new HashMap<String, MSFeature>();
		for(int i = 0; i < this.featureList.size(); i++){
			MSFeature feature1 = featureList.get(i);
			for(int j = i+1; j < this.featureList.size(); j++){
				MSFeature feature2 = featureList.get(j);
				if(checkPair(feature1, feature2, expectedMassDiff, massTolerance, maxRTDiff, mode)){
					System.out.println("Isotopic pairs:\t" + feature1 + "\t" + feature2);
					pairedFeatures.add(new MSFeature[]{feature1, feature2});
					paired.put(feature1.getId(), feature1);
					paired.put(feature2.getId(), feature2);
				}
			}
		}
		
		return pairedFeatures;
	}
	
	public boolean checkPair(MSFeature feature1, MSFeature feature2, double massDiff, double massTolerance, double maxRTDiff, int massMode){
		if(feature1.getCharge() == feature2.getCharge()
				&& Math.abs(feature1.getMeanRT() - feature2.getMeanRT()) < maxRTDiff){
			double mass1 = 0, mass2 = 0;
			if(feature1.getMz() > feature2.getMz()){
				mass2 = feature1.getMz();
				mass1 = feature2.getMz();
			}else{
				mass1 = feature1.getMz();
				mass2 = feature2.getMz();
			}
			mass1 += massDiff/feature1.getCharge();
			return SpectrumUtil.checkMass(mass1, mass2, massTolerance, 2);
		}
		return false;
	}
	
	//get the feature map with time
	public SortedMap<Double, MSFeature> getTimeMap(int mode){
		SortedMap<Double, MSFeature> timeMap = new TreeMap<Double, MSFeature>();
		for(int i = 0; i < featureList.size(); i++){
			MSFeature feature = featureList.get(i);
			if(mode == 0){
				timeMap.put(feature.getMinRT(), feature);
			}else if(mode == 1){
				timeMap.put(feature.getMeanRT(), feature);
			}else{
				timeMap.put(feature.getMaxRT(), feature);
			}
		}
		return timeMap;
	}
	
	//map MS2 spectrum to feature
	public void mapMS2ToFeature(String spectrumFile){
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		System.out.println("Getting spectrum time map");
		SortedMap<Double, Spectrum> sortedRTMap = (SortedMap)reader.getSpectrumTimeMap();
		System.out.println("start mapping MS2");
		for(int i = 0; i < featureList.size(); i++){
			MSFeature feature = featureList.get(i);
			SortedMap sub = sortedRTMap.subMap(feature.getMinRT()-10, feature.getMaxRT()+10);
			//System.out.println("submap size: " + sub.values().size());
			for(Iterator<Spectrum> it = sub.values().iterator(); it.hasNext();){
				Spectrum s = it.next();
				if(feature.isWithnFeature(s.parentMass, s.charge, feature.getMeanRT(), 0.05)){
					if(verbose) System.out.println("mappped spectrum to feature " +  s.scanNumber + "\t" + feature);
					MSFeatureAcq acqFeature = (MSFeatureAcq)feature;
					acqFeature.getAcqInfo().add(s);
				}
			}
		}
		this.mapMS2Yet = true;
		System.out.println("done mapping ms2");
	}
	
	//map IDs to MSfeatures
	//take search results in simple tabullar format separated by tabs
	//note MS2 spectra first has to mappped to feature first
	public void mapIDToFeature(String resultFile, String spectrumFile, int scanInd, int IDInd, int chargeInd){
		if(!this.mapMS2Yet){
			System.err.println("Warning: MS2 scans has not been mapped to MSFeatures");
			return;
		}
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		List<String> results =Utils.FileIOUtils.createListFromFile(resultFile);
		Map<Integer, String> idMap = new HashMap<Integer, String>();
		for(int i = 0; i < results.size(); i++){
			String[] tokens = results.get(i).split("\\t");
			idMap.put(Integer.parseInt(tokens[scanInd]), tokens[IDInd].substring(2, tokens[IDInd].length()-3));
		}
		Map<Integer, String> mapped = new HashMap<Integer, String>();
		Map<String, String> mappedPep = new HashMap<String, String>();
		for(int i = 0; i < featureList.size(); i++){
			MSFeatureAcq acqFeature = (MSFeatureAcq)featureList.get(i);
			List<Spectrum> specList = acqFeature.getAcqInfo();
			for(int j = 0; j < specList.size(); j++){
				Spectrum s = specList.get(j);
				if(idMap.containsKey(s.scanNumber)){
					s.peptide = idMap.get(s.scanNumber);
					mapped.put(s.scanNumber, " ");
					mappedPep.put(s.peptide, "");
					System.out.println("mapped ID: " + s.scanNumber + "\t" + s.peptide + "\t"  
							+ s.charge + "\t" + s.upperBound +  "\tto feature\t" + acqFeature);	
				}
			}
		}
		for(Iterator<Integer> it = idMap.keySet().iterator(); it.hasNext();){
			int scan = it.next();
			String pep = idMap.get(scan);
			Spectrum s = reader.getSpectrum(scan);
			//if(!mapped.containsKey(scan)){
			if(!mappedPep.containsKey(pep)){
				System.out.println("Not-mapped-ID:\t" + scan + "\t" + idMap.get(scan) + "\t" + s.upperBound);
			}
		}
		
	}
	
	public SortedMap<Double, MSFeature> getFeatureMapByMass(){
		SortedMap<Double, MSFeature> featureMap = new TreeMap<Double, MSFeature>();
		List<MSFeature> featureList = this.getFeatureList();
		for(int i = 0; i < this.getFeatureList().size(); i++){
			MSFeature feature = featureList.get(i);
			System.out.println(feature);
			featureMap.put(feature.getMz()*feature.getCharge() - (feature.getCharge()-1)*Mass.PROTON_MASS, feature);
		}
		return featureMap;
	}
	
	public void mapScanFromRT(String spectrumLibFile){
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		List<Spectrum> MS2 = reader.readAllMS2Spectra();
		Map<Integer, Double> rtMap = reader.getRTScanMapping();
		Map<Double, Integer> rtMap2 = reader.getRTScanMappingReverse();
		for(int j = 0; j < this.featureList.size(); j++){
			MSFeature feature = this.featureList.get(j);
			feature.setScan(rtMap2.get(feature.getRt()));
			System.out.println(this.featureList.get(j) );					
		}
	}
	
	public static void testFeatureParser(String featureFile){
		FeatureXMLParser parser = new FeatureXMLParser(featureFile);
	}	
	
	
	public static void testMapMS2ToFeature(){
		//String peptidePairFile = "..//mixture_linked//3UNE_pairs_mis2.txt";
		String spectrumFile = "..//mixture_linked//ProteinQIDs/FeatureDetection//20080313_CPTAC6_16_6D014.mzXML";
		String featureFile = "..//mixture_linked///ProteinQIDs/FeatureDetection//20080313_CPTAC6_16_6D014_binalrecommdParams.FeatureXML";
		String resultFile = "..//mixture_linked//ProteinQIDs/FeatureDetection//20080313_CPTAC6_16_6D014_IDs.txt";
		FeatureXMLParser parser = new FeatureXMLParser(featureFile);
		//double offset = 12.0759;
		//offset = Math.random()*50;
		//System.out.println("offset: " + offset);
		parser.mapMS2ToFeature(spectrumFile);
		parser.mapIDToFeature(resultFile, spectrumFile, 2, 7, 6);
	}
	
	public static void testMapMSFeatureToPeptide(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\Crosslinked_peps_lib-2_cen.mzXML";
		//CombinatoryPeptides combPeps = new CombinatoryPeptides("[NP][ET][GY]K[SQ][IF][AD]R");
		List<String> peptides = combPeps.generateAllPeptides();
		//peptides = addMod(peptides, 42.010564686);
		peptides = combPeps.cutPeptides(peptides);
		System.out.println("we have peptides: " + peptides.size());
		//peptides = addMod(peptides);
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);
		factory.indexPeptideByParentMass(0.1);
		//factory.reIndexPeptideByParentMass(0.1);
		//factory.insertPTM(Mass.DSSDANGLE_MASS, new char[]{'K'}, 1);
		//factory.insertPTM(Mass.DSSDANGLE_MASS_D12, new char[]{'K'}, 1);
		//factory.insertPTM(0.9847, new char[]{'N', 'Q'}, 2);
		factory.crossLinkAllPeptides(2, 4);
		FeatureXMLParser parser = new FeatureXMLParser("..//mixture_linked//test.xml");
		parser.mapScanFromRT(spectrumLibFile);
		List<Peptide> pepList = factory.getAllPeptide();
		List<MSFeature> featureList = parser.getFeatureList();
		for(int i = 0; i < pepList.size(); i++){
			Peptide current = pepList.get(i);
			Peptide currentD12 = new Peptide(current);
			currentD12.setParentmass(currentD12.getParentmass()+12.0759/currentD12.getCharge());
			for(int j = 0; j < featureList.size(); j++){
				MSFeature feature = featureList.get(j);
				if(feature.isMatchFeature(current, 0.03)){
					System.out.println("Peptide: " + current + "\t" + current.getParentmass() + "\t" + current.getCharge() 
							+ " match to feature " + feature);
				}
				if(feature.isMatchFeature(currentD12, 0.03)){
					System.out.println("Peptide-Heavy: " + current + "\t" + current.getParentmass() + "\t" + current.getCharge() 
							+ " match to feature " + feature);
				}
			}
		}
	}
	
	public static void testMapPeptideToMSFeature(){
		CombinatoryPeptides combPeps = new CombinatoryPeptides("[TN][PG][AY]K[EF][IQ][DS]R");
		List<String> peptides = combPeps.generateAllPeptides();
		//peptides = addMod(peptides, 42.010564686);
		peptides = combPeps.cutPeptides(peptides);
		System.out.println("we have peptides: " + peptides.size());
		CandidateSpectrumLibFactory factory = 	
			CandidateSpectrumLibFactory.createFactoryFromPeptide(peptides);
		factory.setMinCharge(1);
		factory.setMaxCharge(1);
		factory.indexPeptideByParentMass(0.1);
		factory.crossLinkAllPeptides(2, 4);
		//String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\20100215_nLCFTMSMS_CLPL_2pmol_1841_MSonly.mzXML";
		String spectrumLibFile = "..\\mixture_linked\\linked_peptide_library\\Crosslinked_peps_lib-2_cen.mzXML";
		MZXMLReader reader = new MZXMLReader(spectrumLibFile);
		PrecursorMassChecker checker = new PrecursorMassChecker(spectrumLibFile);
		long start = (new GregorianCalendar()).getTimeInMillis();
		FeatureXMLParser parser = new FeatureXMLParser("..//mixture_linked//test.xml");
		List<MSFeature> featureList = parser.getFeatureList();
		parser.mapScanFromRT(spectrumLibFile);
		Map<Integer, Double> rtMap = reader.getRTScanMapping();
		List<Peptide> pepList = factory.getAllPeptide();
		for(int i = 0; i < pepList.size(); i++){
			Peptide current = pepList.get(i);
			Peptide currentD12 = new Peptide(current);
			currentD12.setParentmass(currentD12.getParentmass()+12.0759/currentD12.getCharge());
			//System.out.println("checking : " + current + "\t" + current.getCharge());
			//int matched = checker.matchPeptidePrecursorProfile2(current, 10);
			List[] matched = checker.matchPeptidePrecursorProfilePairDetail(current, 70, (Mass.DEUTERIUM-Mass.PROTON_MASS)*12);
			for(int k = 0; k < matched[2].size(); k++){
				System.out.println("checking1 : " + matched[2].get(k) + "\t" + current + "\t" + current.getParentmass()); 
				for(int j = 0; j < featureList.size(); j++){
					MSFeature feature = featureList.get(j);
					double RT = rtMap.get(matched[2].get(k));
					if(feature.isMatchFeature(current, 0.03, RT)){
						System.out.println(current +  "\t" + current.getParentmass() + "\t" + current.getCharge()
								+ " with isotopic labels match to feature " + feature);
						break;
					}
					if(feature.isMatchFeature(currentD12, 0.03, RT)){
						System.out.println(current +  "\t" + currentD12.getParentmass() + "\t" + currentD12.getCharge()
								+ " with isotopic labels match to Heavy-feature " + feature);
						break;
					}
				}				
			}
		}		
		System.out.println("matching spectra in time: " + (new GregorianCalendar().getTimeInMillis()- start)/1000 + "secs");
	}
	
	public static void testMapLinkedPeptideToFeature(){
		//String peptideFile = "..//mixture_linked//database//lib_disulfide2_ACG_pluslengthvariants.txt";
		String peptideFile = "..//mixture_linked//database/lib_disulfide_Ecolidecoy.txt";
		String featureFile = "..//mixture_linked//openMS//feature_detection\\TOPPAS_out\\005-FeatureFinderCentroided\\";
		FeatureXMLParser parser = new FeatureXMLParser(featureFile, "18286_PL2__750ng_5mM_BS3_TRP_SCX_cleanup_IDA.*");
		List<String> peptides = Utils.FileIOUtils.createListFromFile(peptideFile);
		List<Peptide> pepList = new ArrayList<Peptide>(peptides.size());
		for(int k = 0; k < peptides.size(); k++){
			Peptide p = new Peptide(peptides.get(k), 1);
			pepList.add(p);
		}
		
		SortedMap<Double, MSFeature> featureMap = parser.getFeatureMapByMass();

		CrossLinker linker1 = new CrossLinker(138.06, 
				new int[]{CrossLinker.BUTCTERM}, new int[]{CrossLinker.BUTCTERM},
				new char[]{'K'}, new char[]{'K'});
		CrossLinker linker2 = new CrossLinker(-116.43, 
				new int[]{CrossLinker.BUTCTERM}, new int[]{CrossLinker.BUTCTERM},
				new char[]{'C'}, new char[]{'C'});
		CrossLinker linker3 = new CrossLinker(22.025, 
				new int[]{CrossLinker.BUTCTERM}, new int[]{CrossLinker.BUTCTERM},
				new char[]{'C'}, new char[]{'C'});
		CrossLinker linker4 = new CrossLinker(196.117, 
				new int[]{CrossLinker.BUTCTERM}, new int[]{CrossLinker.BUTCTERM},
				new char[]{'C'}, new char[]{'C'});
		CrossLinker linker5 = new CrossLinker(40.0370, 
				new int[]{CrossLinker.BUTCTERM}, new int[]{CrossLinker.BUTCTERM},
				new char[]{'C'}, new char[]{'C'});
		CrossLinker[] linkers = new CrossLinker[]{linker1, linker2, linker3, linker4, linker5};

		for(int i = 0; i < linkers.length; i++){
			List<LinkedPeptide> linkedPeptides = linkers[i].crossLinkPeptides(pepList, 1, 1, true);
			System.out.println("total linked: " + linkedPeptides.size());
			double ppmTolerance = 15;
			Set<MSFeature> mapped = new HashSet<MSFeature>();
			for(int j = 0; j < linkedPeptides.size(); j++){
				Peptide current = linkedPeptides.get(j);
				double deltaMass = current.getParentmass()*ppmTolerance/1000000;
				//System.out.println("deltaMass: " + deltaMass);
				Map subMap = featureMap.subMap(current.getParentmass() - deltaMass, 
						current.getParentmass() + deltaMass);
				if(subMap.size() > 0){
					System.out.println("macthed size: " + subMap.size());
					System.out.println("mapping peptide " + current + " to feature " + subMap.values().iterator().next());
				}
				mapped.addAll(subMap.values());
			}
			System.out.println("Linked peptide " + i + " mapped to features: " + mapped.size());
		}
	}
	
	
	public static void testMapPeptidePairToMSFeature(){
		String peptidePairFile = "..//mixture_linked//3UNE_pairs_mis2.txt";
		String spectrumFile = "..//mixture_linked//msdata//proteasome_crosslinks//aleitner_M1108_136.mzXML";
		String featureFile = "..//mixture_linked//openMS//feature_detection\\TOPPAS_out\\005-FeatureFinderCentroided\\";
		String resultFile = "..//mixture_linked//testAnnotation.txt1";
		List<String> peptides = FileIOUtils.createListFromFile(peptidePairFile);
		List<Peptide> pepList = new ArrayList<Peptide>();
		FeatureXMLParser parser = new FeatureXMLParser(featureFile, "aleitner_M1108_136.*XML");
		SortedMap<Double, MSFeature> featureMap = parser.getFeatureMapByMass();
		Map<Peptide, String> siteMap = new HashMap<Peptide, String>();
		parser.mapMS2ToFeature(spectrumFile);
		parser.mapIDToFeature(resultFile, spectrumFile, 1, 7, 6);
		for(int i = 0; i < peptides.size(); i++){
			String[] tokens = peptides.get(i).split("\\s+");
			double linkerMass = 138.068;
			//linkerMass = Math.random()*13;
			LinkedPeptide lp = new LinkedPeptide(new Peptide(tokens[2],1), new Peptide(tokens[3],1), 1, 
					tokens[2].indexOf('K')+1, tokens[3].indexOf('K')+1, linkerMass);
			//Peptide lp = new Peptide(tokens[2],1);
			pepList.add(lp);
			siteMap.put(lp, tokens[5]);
		}
		
		System.out.println("we have peptides: " + peptides.size());
		double tolerance = 15;
		for(int i = 0; i < pepList.size(); i++){
			Peptide current = pepList.get(i);
			Peptide currentD12 = new Peptide(current);
			String site = siteMap.get(current);
			System.out.println("peptide: " + current);
			currentD12.setParentmass(currentD12.getParentmass()+12.0759/currentD12.getCharge());
			double deltaMass = current.getParentmass()*tolerance/1000000;
			System.out.println("deltaMass: " + deltaMass);
			Map subMap = featureMap.subMap(current.getParentmass() - deltaMass, 
					current.getParentmass() + deltaMass);
			if(subMap.size() > 0){
				System.out.println("macthed size: " + subMap.size());
				System.out.println("mapping peptide " + current + " to feature " + subMap.values().iterator().next() +"\t" + site);
			}
			Map subMap2 = featureMap.subMap(currentD12.getParentmass() - deltaMass, 
					currentD12.getParentmass() + deltaMass);
			if(subMap2.size() > 0){
				System.out.println("macthed size: " + subMap.size());
				System.out.println("mapping peptide D12 " + current + " to feature " + subMap2.values().iterator().next() +"\t" + site);
			}
			if(subMap.size() > 0 && subMap2.size() > 0){
				for(Iterator it1 = subMap.values().iterator(); it1.hasNext();){
					MSFeatureAcq feature = (MSFeatureAcq)it1.next();
					for(Iterator it2 = subMap2.values().iterator(); it2.hasNext();){
						MSFeatureAcq feature2 = (MSFeatureAcq)it2.next();
						if(parser.checkPair(feature, feature2, 12.0759, 30, 30, 2)){
							System.out.println("mapping peptide " + current + " to feature-pair " + feature + " & " + feature2 +"\t" + site 
									+"\tacquired:\t" + feature.getAcqInfo().size() + "\t" + feature2.getAcqInfo().size());
						}
					}
				}
			}
			//System.out.println("checking : " + current + "\t" + current.getCharge());
			//int matched = checker.matchPeptidePrecursorProfile2(current, 10);
			
		}	
	}
	
	public static void testFindFeaturePair(){
		String peptidePairFile = "..//mixture_linked//3UNE_pairs_mis2.txt";
		String spectrumFile = "..//mixture_linked//msdata//proteasome_crosslinks//aleitner_M1108_136.mzXML";
		String featureFile = "..//mixture_linked//openMS//feature_detection\\TOPPAS_out\\005-FeatureFinderCentroided\\";
		String resultFile = "..//mixture_linked//testAnnotation.txt1";
		FeatureXMLParser parser = new FeatureXMLParser(featureFile, "aleitner_M1108_136.*");
		double offset = 12.0759;
		//offset = Math.random()*50;
		System.out.println("offset: " + offset);
		parser.mapMS2ToFeature(spectrumFile);
		parser.mapIDToFeature(resultFile, spectrumFile, 2, 7, 6);
		List<MSFeature[]> pairs = parser.getFeaturePair(offset, 30, 30, 2);
		for(int i = 0; i < pairs.size(); i++){
			MSFeature[] featurePair = pairs.get(i);
			if(featurePair[0].getId().equals("f_10738625050898221243") || featurePair[1].getId().equals("f_10738625050898221243")){
				//System.out.println("here");
			}
			List<Spectrum> msms1 = ((MSFeatureAcq)featurePair[0]).getAcqInfo();
			for(int j = 0; j < msms1.size(); j++){
				Spectrum s = msms1.get(j);
				System.out.println("pair: " + featurePair[0] + "\tmsms:\t" +  s.scanNumber  + "\t" + s.peptide +"\t" + s.charge);
			}
			List<Spectrum> msms2 = ((MSFeatureAcq)featurePair[1]).getAcqInfo();
			//System.out.println("acquired list: " + msms1.size() + "\t" + msms2.size());
			for(int j = 0; j < msms2.size(); j++){
				Spectrum s = msms2.get(j);
				System.out.println("pair: " + featurePair[1] + "\tmsms:\t" +  s.scanNumber  + "\t" + s.peptide + "\t" + s.charge);
			}
		}
	}
	
	
	
	public static void main(String[] args){
		String featureFile = "..\\mixture_linked\\ProteinQIDs/FeatureDetection\\CPTACT_featureFinder_0.featureXML";
		//String featureFile = args[0];
		//testFeatureParser(featureFile);
		testMapMS2ToFeature();
		//testMapMSFeatureToPeptide();
		//testMapPeptideToMSFeature();
		//testFindFeaturePair();
		//testMapPeptidePairToMSFeature();
		// testMapLinkedPeptideToFeature();
	}
}

