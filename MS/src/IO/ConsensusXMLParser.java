package IO;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.Spectrums.MSFeature;
import org.Spectrums.MSFeatureAcq;
import org.Spectrums.ProteinIDExtractor;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import Utils.ArrayUtils;

/**
 * Parse consensusXML from openMS
 * @author Jian
 *
 */
public class ConsensusXMLParser {
	
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
	public ConsensusXMLParser(String filename){
		File path = new File(filename);
//		if(path.isDirectory()){
//			File[] files = path.listFiles();
//			for(int i = 0; i < files.length; i++){
//				String filepath = files[i].getAbsolutePath();
//				FeatureXMLParser parser = new FeatureXMLParser(filepath);
//				this.featureList.addAll(parser.getFeatureList());
//			}
//		}else if(path.exists()){
			parseXmlFile(filename);
			parseDocument();
//		}
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
		String fasta = "..//mixture_linked//database//UPS_plusHuman_plusDecoy.fasta";
		ProteinIDExtractor protID = getProteins(docEle, fasta);
		Map<String, List<String>> pepMap = protID.getPeptideMap();
		NodeList nl = docEle.getElementsByTagName("consensusElement");
		System.out.println("There are total of features: " + nl.getLength());
//		this.featureList = new ArrayList();
		//ProteinIDExtractor protID = new ProteinIDExtractor()
		System.out.println("#FeatureID\t#-aligned\tAvgIt\tItVar\tItCV\t#IDs-mapped\tConsistency");
		if(nl != null && nl.getLength() > 0) {
			int count = 0;
			for(int i = 0 ; i < nl.getLength(); i++) {
				Element el = (Element)nl.item(i);
				NodeList elements = el.getElementsByTagName("element");
				NodeList ids = el.getElementsByTagName("PeptideHit");
				List<Double> abundances = new ArrayList<Double>();
				List<String> peptides = new ArrayList<String>();
				for(int j = 0; j < elements.getLength(); j++){
					Element el2 = (Element)elements.item(j);
					//System.out.println("Total intensity: " + el2.getAttribute("it"));
					abundances.add(Double.parseDouble(el2.getAttribute("it")));
				}
				for(int k = 0; k < ids.getLength(); k++){
					Element el3 = (Element)ids.item(k);
					//System.out.println("Total intensity: " + el2.getAttribute("it"));
					//System.out.println("ID: " + el3.getAttribute("sequence"));
					peptides.add(el3.getAttribute("sequence"));
				}

				double consistency = checkPeptides(peptides);
				peptides.add("PEPTIDE");
				double[] stat = getStat(abundances);
				List<String> protList = pepMap.get(peptides.get(0));
				if(protList == null){
					protList=new ArrayList<String>();
					protList.add("PROTEIN");
				}
				System.out.println("Feature: " + el.getAttribute("id") + "\t" + abundances.size() + "\t" + stat[0] + "\t" + stat[1] + "\t" + stat[2] + "\t" 
						+peptides.size() + "\t" + consistency + "\t" + peptides.get(0) + "\t" + protList.get(0) + "\t" + protList.size() +"\t" + el.getAttribute("charge"));
				
			}
		}
	}
	
	public ProteinIDExtractor getProteins(Element docEle, String fastaFile){
		NodeList ids = docEle.getElementsByTagName("PeptideHit");
		Set<String> peptides = new HashSet<String>();
		for(int i = 0; i < ids.getLength(); i++){
			Element el3 = (Element)ids.item(i);
			peptides.add(el3.getAttribute("sequence"));
		}
		System.out.println("peptide length " + peptides.size());
		ProteinIDExtractor protID = new ProteinIDExtractor(peptides, fastaFile);
		return protID;
	}
	
	public double checkPeptides(List<String> peptides){
		int matches = 0;
		int total = 0;
		for(int i = 0; i < peptides.size(); i++){
			String pep1 = peptides.get(i);
			for(int j = i+1; j < peptides.size(); j++){
				String pep2 = peptides.get(j);
				if(pep1.equals(pep2)){
					matches++;
				}
				total++;
			}
		}
		if(total == 0){
			return 1;
		}
		//System.out.println("IDs consisitency: " + matches/total);
		return matches/total;
	}
	
	public double[] getStat(List<Double> values){
		//double mean = ArrayUtils.average(values.to)
		double[] arry1 = ArrayUtils.toArray(values);
		double avg = ArrayUtils.average(arry1);
		double var = ArrayUtils.var(arry1);
		double CV = ArrayUtils.CV(arry1);
		return new double[]{avg, var, CV};
	}
	
	public static void testParseConsensusDocument(){
		String file = "..//mixture_linked/ProteinQIDs/OpenMS/TOPPAS_out/TOPPAS_out/020-FeatureLinkerUnlabeledQT-out/database_tmp25.consensusXML";
		//String file = "..//mixture_linked/ProteinQIDs//OpenMS/CPTAC_FeatureLinkerUnlabeledQT_technical_rep/yeast_proteins_plusDecoy.consensusXML";
		//String file = "..//mixture_linked/ProteinQIDs//OpenMS/CPTAC_FeatureLinkerUnlabeledQT_bio_rep/017-FeatureLinkerUnlabeledQT/yeast_proteins_plusDecoy.consensusXML";
		ConsensusXMLParser parser = new ConsensusXMLParser(file);
	}
	
	public static void main(String[] args){
		testParseConsensusDocument();
	}

}
