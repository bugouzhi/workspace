import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Large hashmap implementation back-by a random access file
 * @author Jian Wang
 *
 */
public class LargeHashMap {
	private String libraryObjectFile;
	private int keepObjects = 20; //number of objects keep in memory
	public boolean DEBUG=false;
	public int getKeepObjects() {
		return keepObjects;
	}


	public void setKeepObjects(int keepObjects) {
		this.keepObjects = keepObjects;
	}

	private RandomAccessFile RAF;
	private Map table;
	private SortedSet keys;
	private Map<Object, Long> beginTable;
	private Map<Object, Long> endTable;
	
	//load from table files
	public LargeHashMap(String libraryObjectFile){
		this.libraryObjectFile = libraryObjectFile;
		this.table = new HashMap();
		this.keys = new TreeSet();
		this.loadLibraryFromFile(this.libraryObjectFile);
	}
	
	//build the table
	public LargeHashMap(Map table, String libraryObjectFile){
		this.libraryObjectFile = libraryObjectFile;
		this.table = new HashMap();
		this.keys = new TreeSet();
		buildTable(table);
		this.loadLibraryFromFile(this.libraryObjectFile);
	}
	
	
	public void loadLibraryFromFile(String file){
		this.libraryObjectFile = file;
		try{
			this.RAF = new RandomAccessFile(file, "r");
			long index1 = RAF.readLong();
			long index2 = RAF.readLong();
			long index3 = RAF.readLong();
			long current = RAF.getFilePointer();
			if(DEBUG)
				System.out.println("index is: " + index1 + "\t" + index2 + "\t" + index3);
			RAF.seek(index1);
			byte[] input = new byte[(int)(index2-index1)];
			System.out.println("current is : " + RAF.getFilePointer());
			int readed = RAF.read(input);
			if(DEBUG)	
				System.out.println("read in bytes: " + readed);
			ByteArrayInputStream in = new ByteArrayInputStream(input);
			ObjectInputStream in2 = new ObjectInputStream(in);
			this.beginTable = (Map)in2.readObject();
			in.close();
			in2.close();
			RAF.seek(index2);
			input = new byte[(int)(index3-index2)];
			if(DEBUG)
				System.out.println("current is : " + RAF.getFilePointer());
			readed = RAF.read(input);
			if(DEBUG)
				System.out.println("read in bytes: " + readed);
			in = new ByteArrayInputStream(input);
			in2 = new ObjectInputStream(in);
			this.endTable = (Map)in2.readObject();
			in.close();
			in2.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}catch(ClassNotFoundException cne){
			System.out.println(cne.getMessage());
			cne.printStackTrace();
		}

	}
	
	public void buildTable(Map specTable){
		try{
			this.beginTable = new HashMap<Object, Long>();
			this.endTable = new HashMap<Object, Long>();
			RandomAccessFile RAF = new RandomAccessFile(this.libraryObjectFile, "rw");
			//first we store file position to beginTable and endTable
			long index1 = RAF.getFilePointer();
			RAF.writeLong(1L);  
			long index2 = RAF.getFilePointer();
			RAF.writeLong(1L);
			long index3 = RAF.getFilePointer();
			RAF.writeLong(1L);  
			
			ByteArrayOutputStream out1;
			ObjectOutputStream out2;
			
			for(Iterator it = specTable.keySet().iterator(); it.hasNext();){
				int key = (Integer)it.next();
				beginTable.put(key, RAF.getFilePointer());
				//System.out.println("begin is: " + RAF.getFilePointer());
				Object values = specTable.get(key);
				out1 = new ByteArrayOutputStream();
				out2 = new ObjectOutputStream(out1);
				out2.writeObject(values);
				out2.flush();
				RAF.write(out1.toByteArray());
				endTable.put(key, RAF.getFilePointer());
				//System.out.println("end is : " + RAF.getFilePointer());
			}
			
			//storing the index table
			long current = RAF.getFilePointer();
			RAF.seek(index1);
			RAF.writeLong(current);
			RAF.seek(current);
			if(DEBUG)
				System.out.println("begin writing " + current);
			out1 = new ByteArrayOutputStream();
			out2 = new ObjectOutputStream(out1);
			out2.writeObject(this.beginTable);
			out2.flush();
			RAF.write(out1.toByteArray());
			if(DEBUG)
				System.out.println("end writing " + RAF.getFilePointer() + " supposely written: " + out1.toByteArray().length);
			out1.close();
			out2.close();
			
			current = RAF.getFilePointer();
			RAF.seek(index2);
			RAF.writeLong(current);
			RAF.seek(current);
			
			if(DEBUG)
				System.out.println("begin writing " + current);
			out1 = new ByteArrayOutputStream();
			out2 = new ObjectOutputStream(out1);
			out2.writeObject(this.endTable);
			out2.flush();
			RAF.write(out1.toByteArray());
			if(DEBUG)
				System.out.println("end writing " + RAF.getFilePointer() + " supposely written: " + out1.toByteArray().length);
			current = RAF.getFilePointer();
			RAF.seek(index3);
			RAF.writeLong(current);
			RAF.seek(current);
			out1.close();
			out2.close();
			//out2.writeObject(this.beginTable);
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public Object get(Object key){
		if(table.containsKey(key)){
			return this.table.get(key);
		}else if(this.beginTable.containsKey(key)){
			try{
				long begin = this.beginTable.get(key);
				long end = this.endTable.get(key);
				RAF.seek(begin);
				byte[] input = new byte[(int)(end-begin)];
				RAF.read(input);
				ByteArrayInputStream in = new ByteArrayInputStream(input);
				ObjectInputStream in2 = new ObjectInputStream(in);
				Object value = in2.readObject();
				//System.out.println("Object is: " + value.getClass());
				if(this.keys.size() == this.keepObjects){
					Object smallest = keys.first();
					this.keys.remove(smallest);
					this.table.remove(smallest);
				}
				this.table.put(key, value);
				this.keys.add(key);
				return value;
			}catch(IOException ioe){
				System.out.println(ioe.getMessage());
				ioe.printStackTrace();
			}catch(ClassNotFoundException cne){
				System.out.println(cne.getMessage());
				cne.printStackTrace();
			}
			return null;
		}else{
			return null;
		}
	}

}
