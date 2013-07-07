package Utils;
/**
 * Contain utility function for performing various file-related i/o task
 * @author jian wang
 *
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Iterator;

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
		List<String> lines = new ArrayList<String>();
		try{
			BufferedReader bf = new BufferedReader(new FileReader(file));
			String line = bf.readLine();
			//lines.add(line);
			while(line != null){
				lines.add(line);
				line = bf.readLine();
			}
			bf.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
		}catch(NullPointerException e){
			System.out.println("Error reading in file: " + file);
			System.out.println(e.getMessage());
			System.out.println(e.getStackTrace());
		}
		System.out.println("readed in " + lines.size() + " lines");
		return lines;
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

	public static String stripExtension(String filename){
		int index = filename.lastIndexOf(".");
		return filename.substring(0, index);
	}
	
}
