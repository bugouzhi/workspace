package UI;
import org.Spectrums.ResultColumnIndex;
import org.Spectrums.SpectrumLibConstructor;
import org.Spectrums.SpectrumUtil;

import Utils.SpectrumFilter;
/**
 * Create a annoated mgf from raw file and search results
 * This is a wrapper class for running a specific method from the command line
 * @author Jian
 *
 */
public class CreateSpectralLibrary {
	public static void testCreateSpecLib(){
		String spectrumDir = "../mixture_linked/Emily_toni_Exosomes/mgfs/";
		String annotationFile = "../mixture_linked/Emily_toni_Exosomes/Mouse_exosome_DDA_mgfs_msgfdb_results_sorted.txt";
		String outfile = "..\\mixture_linked\\test.mgf";
		String logfile = "../mixture_linked/creatLib.log";
		runCreateSpectLib(new String[]{spectrumDir, annotationFile, outfile, logfile});
	}
	
	//legacy method using spectrumUtil's annotateSpectrumLibFromMzXML rather than the spectrumlibConstructor class
	public static void runCreateSpectLib_le(String[] args){
		CommandLineParser parser = new CommandLineParser(args);
		if(args.length == 4){
			SpectrumUtil.annotateSpectrumLibFromMzXMLs(args[0], args[1], args[2], args[3]);
		}else if(args.length == 10){
			int mode = parser.getInteger(4);
			SpectrumFilter filter = new SpectrumFilter();
			if(mode == SpectrumFilter.MODE_GLOBAL){
				double lowFract = parser.getDouble(5);
				filter.setFilteringMode(mode);
				filter.setLowFract(lowFract);
			}else{
				filter.setTopN(parser.getInteger(6));
				filter.setWindowWidth(parser.getDouble(7));
			}
			filter.setFilteringMode(mode);
			filter.setMinSignalR(parser.getDouble(8));
			double tolerance = parser.getDouble(9);
			SpectrumUtil.annotateSpectrumLibFromMzXMLs(args[0], args[1], args[2], args[3], filter, tolerance);
		}else{
			System.out.println("java -cp MSPLIT-DIA.jar UI.CreateSpectralLibrary <SpectrumDir> <annotationFile> <mgf out file> <log file>\n" + 
					"optional <filtering mode (0: global | 1:local)> <fraction least-intense peaks> <topN> <window width> <min SignalRatio> <fragment mass tolerance>");
		}
	}
	
	public static void runCreateSpectLib(String[] args){
		CommandLineParser parser = new CommandLineParser(args);
		SpectrumLibConstructor constructor = new SpectrumLibConstructor(parser.getString(0), parser.getString(1), ResultColumnIndex.MSGFDB_INDEX);
		if(args.length == 4){
			SpectrumFilter filter = new SpectrumFilter();
			constructor.constructSpectralLibrary(parser.getString(2), parser.getString(2), filter, 0.05);
		}else if(args.length == 10){
			int mode = parser.getInteger(4);
			SpectrumFilter filter = new SpectrumFilter();
			double lowFract = parser.getDouble(5);
			filter.setLowFract(lowFract);
			filter.setTopN(parser.getInteger(6));
			filter.setWindowWidth(parser.getDouble(7));
			filter.setMinSignalR(parser.getDouble(8));
			if(mode == 0){
				filter.setGlobalFilter(true);
			}else if(mode == 1){
				filter.setLocalFilter(true);
			}else if(mode == 2){
				filter.setGlobalFilter(true);
				filter.setLocalFilter(true);
			}
			double tolerance = parser.getDouble(9);
			constructor.constructSpectralLibrary(parser.getString(2), parser.getString(3), filter, tolerance);
		}else{
			System.out.println("java -cp MSPLIT-DIA.jar UI.CreateSpectralLibrary <SpectrumDir> <annotationFile> <mgf out file> <log file>\n" + 
					"optional <filtering mode (0: global | 1:local | 2: both)> <fraction least-intense peaks> <topN> <window width> <min SignalRatio> <fragment mass tolerance>");
		}
	}
	
	public static void main(String[] args){
		//testCreateSpecLib();
		runCreateSpectLib(args);
	}
}




