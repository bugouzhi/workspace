package Utils;
/**
 * Contain utility function for performing various file-related i/o task
 * @author jian wang
 *
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Iterator;
//import org.biojava.bio.Annotation;
//import org.biojava.bio.BioException;
//import org.biojava.bio.seq.ProteinTools;
//import org.biojava.bio.seq.io.CharacterTokenization;
//import org.biojava.bio.seq.io.SymbolTokenization;
//import org.biojava.bio.symbol.AlphabetManager;
//import org.biojava.bio.symbol.FiniteAlphabet;
//import org.biojava.bio.symbol.SimpleAlphabet;
//import org.biojava.bio.symbol.Symbol;
//import org.biojava.utils.ChangeVetoException;
//import org.biojavax.RichObjectFactory;
//import org.biojavax.bio.db.HashRichSequenceDB;
//import org.biojavax.bio.seq.RichSequence;
//import org.biojavax.bio.seq.RichSequenceIterator;
//import org.biojavax.bio.seq.io.RichSequenceBuilderFactory;

public class FileIOUtils {
	public static Map<String, String> createTableFromFile(String file, int index1, int index2){
		Map<String, String> table = new LinkedHashMap<String, String>();
		if(index1 < 0 || index2 < 0){
			throw new IllegalArgumentException("Index is less than zero");
		}
		int max = index1 > index2 ? index1 : index2;
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			String[] token = null;
			String line = bf.readLine();
			while(line != null){
				System.out.println("line is: " + line);
				token = line.split("\\t");
				table.put(token[index1], token[index2]);
				line = bf.readLine();
			}
			bf.close();
			System.out.println("reading in " + table.keySet().size() + " lines from file");
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
		}catch(NullPointerException e){
			System.out.println("Error reading in file: " + file);
			System.out.println(e.getMessage());
			System.out.println(e.getStackTrace());
		}
		return table;
	}
	
	public static List<String> createListFromFile(String file){
		return readLines(file, Integer.MAX_VALUE);
	}
	
	
	public static List<String> createListFromFile(String file, int index, String separator){
		return readLines(file, Integer.MAX_VALUE, index, separator);
	}
	
	public static List<String> readLines(String file, int num, int index, String separator){
		List<String> lines = new ArrayList<String>();
		int count = 0;
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			String line = bf.readLine();
			//lines.add(line);
			while(line != null && count < num){
				String[] tokens = line.split(separator);
				lines.add(tokens[index]);
				line = bf.readLine();
				count++;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
		}catch(NullPointerException e){
			System.out.println("Error reading in file: " + file);
			System.out.println(e.getMessage());
			System.out.println(e.getStackTrace());
		}
		//System.out.println("readed in " + lines.size() + " lines");
		return lines;
	}
	
	public static List<String> readLines(String file, int num){
		List<String> lines = new ArrayList<String>();
		int count = 0;
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			String line = bf.readLine();
			//lines.add(line);
			while(line != null && count < num){
				lines.add(line);
				line = bf.readLine();
				count++;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
		}catch(NullPointerException e){
			System.out.println("Error reading in file: " + file);
			System.out.println(e.getMessage());
			System.out.println(e.getStackTrace());
		}
		//System.out.println("readed in " + lines.size() + " lines");
		return lines;
	}
	
	public static List<String> readLines(String file, int num, boolean check){
		List<String> lines = new ArrayList<String>();
		int count = 0;
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			String line = bf.readLine();
			//lines.add(line);
			while(line != null && count < num){
				lines.add(line);
				line = bf.readLine();
				count++;
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
		}catch(NullPointerException e){
			System.out.println("Error reading in file: " + file);
			System.out.println(e.getMessage());
			System.out.println(e.getStackTrace());
		}
		//System.out.println("readed in " + lines.size() + " lines");
		return lines;
	}
	
	public static BufferedReader createReaderFromFile(String file){
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			return bf;
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
		}catch(NullPointerException e){
			System.out.println("Error reading in file: " + file);
			System.out.println(e.getMessage());
			System.out.println(e.getStackTrace());
		}
		return null;
		//System.out.println("readed in " + lines.size() + " lines");
	}
	
	
	public static List<String> createProteinsFromFasta(String file){
		List<String> lines = createListFromFile(file);
		List<String> proteins = new ArrayList<String>();
		StringBuffer buffer = new StringBuffer();
		for(int i = 0; i < lines.size(); i++){
			String line = lines.get(i);
			if(line.startsWith(">")){
				if(buffer.length() > 0){
					proteins.add(buffer.toString());
					buffer = new StringBuffer();
				}
			}else{
				buffer.append(line);
			}
		}
		return proteins;
	}
	
	
	public static String stripExtension(String filename){
		int index = filename.lastIndexOf(".");
		return filename.substring(0, index);
	}
	
	
	public static String getFileExtension(String filename){
		int index = filename.lastIndexOf(".");
		if(index > 0 && index < filename.length()){
			return filename.substring(index+1);
		}else{
			return "";    //return empty string if there is no extension
		}
	}
	
	//get only the file name from the path 
	public static String getFileName(String path){
		int ind = -1;
		if(path.contains("/")){
			ind = path.lastIndexOf("/");
		}
		if(path.contains("\\")){
			ind = path.lastIndexOf("\\");
		}
		//System.out.println("Filename is: " + path.substring(ind+1));
		return path.substring(ind+1);
	}
	
	public static void createFileFromList(List<String> lines, String file){
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			for(int i = 0; i < lines.size(); i++){
				String line = lines.get(i);
				bw.append(line);
				bw.append("\n");
				//lines.add(line);
			}
			bw.flush();
			bw.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
		}catch(NullPointerException e){
			System.out.println("Error writing to file: " + file);
			System.out.println(e.getMessage());
			System.out.println(e.getStackTrace());
		}
		//System.out.println(" in " + lines.size() + " lines");
	}
	
	//a convinient method to print output, if specify a valid file
	//print to a outfile otherwise print to stdout
	public static PrintStream getOutStream(String outfile){
		try{
			File f = new File(outfile);
			PrintStream out = new PrintStream(f);
			return out;
		}catch(IOException ioe){
			System.out.println("Using standard out as output");
			System.out.println(ioe.getMessage());
			return System.out;
		}
	}
	
	public static void mergeResultWithHeader(String dir, String merged){
		File resultDir = new File(dir);
		if(!resultDir.exists()){
			System.err.println(dir + " do not exits");
			return;
		}
		int counter = 0;
		File[] files = resultDir.listFiles();
		System.out.println("Total number of files: " + files.length);
		PrintStream out = Utils.FileIOUtils.getOutStream(merged);
		for(int i = 0; i < files.length; i++){
			List<String> lines = Utils.FileIOUtils.createListFromFile(files[i].getAbsolutePath());
			for(int j = 0; j < lines.size(); j++){
				String line = lines.get(j);
				if(!(line.startsWith("#") && j == 0)){
					out.println(line);
					counter++;
				}
			}
		}
		System.out.println("Combining files, total lines: " + counter);
		out.flush();
		out.close();
		
		
	}
	
	/**
	 * check if a path is a directory
	 * @param filename
	 * @return
	 */
	public static boolean checkDirectory(String filename){
		File f = new File(filename);
		return f.isDirectory();
	}
	
	/**
	 * Generate a new file name base on another file name
	 * by stripping the extension and adding new extension
	 */
	public static String rename(String oldPathName, String extension){
		String name = stripExtension(oldPathName);
		return name+extension;
	}
	
	/**
	 * Generate output file based on the given output path, if it is a valid file
	 * use that file, if it is a directory, then generate a appropriate file name base on 
	 * the name of another file by stripping the extension of old file and adding another extension
	 * @param out
	 * @param oldName
	 * @param extension
	 * @return
	 */
	public static String generateOutFile(String out, String oldName, String extension){
		if(checkDirectory(out)){
			return out+File.separator+ stripExtension(getFileName(oldName)) + extension;
		}else{
			return out;
		}
	}
	
	
//	public static HashRichSequenceDB  readAlnFile(String file){
//        HashRichSequenceDB db = new HashRichSequenceDB();
//		try{
//			BufferedReader bf = new BufferedReader(new FileReader(file));
//			FiniteAlphabet extendedAlphabet = createCustomAlphabet(); //standard protein alphabet
//			SymbolTokenization rParser = extendedAlphabet.getTokenization("token");
//			
//			RichSequenceIterator seqI = RichSequence.IOTools.readFasta(bf, rParser, RichSequenceBuilderFactory.THRESHOLD, 
//					RichObjectFactory.getDefaultNamespace());
//			int counter = 0;
//			System.out.println("done reading file");
//			while (seqI.hasNext()) {
//                try {
//                    db.addRichSequence(seqI.nextRichSequence());
//                  // System.out.println("readed in " + ++counter + " sequences");
//                    
//                } catch (ChangeVetoException ce) {
//                    throw new BioException("Unexpectedly couldn't add to the supplied RichSequenceDB", ce);
//                }
//            }
//		}catch(Throwable t){
//			System.err.println("unable to read the multiple seq alignments");
//			t.printStackTrace();
//		}
//		return db;
//	}
//	
//	private static FiniteAlphabet createCustomAlphabet(){
//		 	Symbol maskx = AlphabetManager.createSymbol("maskx", Annotation.EMPTY_ANNOTATION);
//			FiniteAlphabet proteinAlphabet = ProteinTools.getTAlphabet(); //standard protein alphabet
//			SimpleAlphabet extendedProteinAlphabet = new SimpleAlphabet();
//			try{
//				for(Iterator<Symbol> it = proteinAlphabet.iterator(); it.hasNext();){
//					extendedProteinAlphabet.addSymbol(it.next());
//				}
//				extendedProteinAlphabet.addSymbol(maskx);
//				SymbolTokenization rParser = proteinAlphabet.getTokenization("token");
//				CharacterTokenization eParser = new CharacterTokenization(extendedProteinAlphabet, true);
//				for(Iterator<Symbol> it = proteinAlphabet.iterator(); it.hasNext();){
//					Symbol c = it.next();
//					eParser.bindSymbol(c, rParser.tokenizeSymbol(c).charAt(0));
//				}
//				eParser.bindSymbol(maskx, 'x');
//				extendedProteinAlphabet.putTokenization("token", eParser);
//				
//			}catch(Exception ioe){
//				
//			}
//			
//			return extendedProteinAlphabet;
//	}
	
	
	
}
