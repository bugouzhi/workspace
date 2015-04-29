package Utils;

import java.io.File;
import java.io.FilenameFilter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Combine multiple tables into one.
 * User specify which column in the table is use as key.
 * @author Jian
 *
 */
public class TableCombiner {
	private String dir;
	private String missingValue = "0";
	private PrintStream out = System.out;
	public TableCombiner(String dir){
		this.dir = dir;
	}
	
	public void combineTable(int keyInd, int valueInd, String outFile){
		this.out = FileIOUtils.getOutStream(outFile);
		File directory = new File(dir);
		String[] files = directory.list(new FilenameFilter(){
			@Override
			public boolean accept(File arg0, String arg1) {
				return arg1.endsWith(".txt");
			}
			
		});
		List<Map<String, String>> tables = new ArrayList<Map<String, String>>();
		Set<String> masterKeySet = new HashSet<String>();
		for(int i = 0; i < files.length; i++){
			String path = this.dir + File.separator + files[i];
			Map<String, String> table = FileIOUtils.createTableFromFile(path, keyInd, valueInd);
			masterKeySet.addAll(table.keySet());
			tables.add(table);
		}
		//printing header
		out.print("#Key\t");
		for(int i = 0; i < tables.size(); i++){
			String name = files[i].split("_")[2];
			out.print(name + "\t");
		}
		out.println();
		for(Iterator<String> it = masterKeySet.iterator(); it.hasNext();){
			String key = it.next();
			out.print(key + "\t");
			for(int i = 0; i < tables.size(); i++){
				Map<String,String> table = tables.get(i);
				if(table.containsKey(key)){
					out.print(table.get(key));
				}else{
					out.print(this.missingValue);
				}
				out.print("\t");
			}
			out.println();
		}
	}
	
	public static void testCombineTable(){
		String dir = "..//mixture_linked//Emily_Toni_DIA_vs_DDA/WiSIM/SpecCount";
		TableCombiner combine = new TableCombiner(dir);
		combine.combineTable(2, 3, "..//mixture_linked//tableCombine.txt");
	}
	
	public static void main(String[] args){
		testCombineTable();
	}

}
