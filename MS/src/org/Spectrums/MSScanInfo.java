package org.Spectrums;

import java.util.Iterator;
import java.util.List;

import org.systemsbiology.jrap.stax.Scan;

import IO.MZXMLReader;
import Utils.FileIOUtils;

/**
 * Add MS scan level information to IDs
 * @author Jian Wang
 *
 */
public class MSScanInfo {
	public static void addMSLevelInfo(){
		String spectrumFile = "../mixture_linked/msdata/philAndrews/Aldolase_BS3_032408/Aldolase_xlink_BS3_032508_scx7.mzXML";
		String resultFile = "../mixture_linked/MSGFDBSearches/aldolase_scx7_mzxml.out";
		List<String> results = FileIOUtils.createListFromFile(resultFile);
		MZXMLReader reader = new MZXMLReader(spectrumFile);
		for(Iterator<String> it = results.iterator(); it.hasNext();){
			String line = it.next();
			if(line.contains("#")){
				continue;
			}
			String[] tokens = line.split("\\t");
			int scan = Integer.parseInt(tokens[1]);
			Scan s = reader.parser.rap(3009);
			//System.out.println(s.getHeader().getMsLevel() +"\t" + s.getHeader().getPrecursorMz());
			System.out.println(line + "\t" + reader.parser.rap(scan).getHeader().getMsLevel());
		}
	}
	
	public static void main(String[] args){
		addMSLevelInfo();	
	}
}
