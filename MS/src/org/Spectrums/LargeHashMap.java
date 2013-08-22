package org.Spectrums;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
	private int keepObjects = 30; //number of objects keep in memory
	private RandomAccessFile RAF;
	private Map table;
	private SortedSet keys;
	private Map<Object, Long> beginTable = new HashMap();
	private Map<Object, Long> endTable = new HashMap();
	private long index1=-1;
	private long index2=-1;
	private long index3=-1;
	private long lastEntryIndex=0;
	public boolean DEBUG=false;
	public int getKeepObjects() {
		return keepObjects;
	}


	public void setKeepObjects(int keepObjects) {
		this.keepObjects = keepObjects;
	}

	
	//build the table
	public LargeHashMap(Map table, String libraryObjectFile){
		this.libraryObjectFile = libraryObjectFile;
		this.table = new HashMap();
		this.keys = new TreeSet();
		buildTable(table);
		this.loadLibraryFromFile(this.libraryObjectFile);
	}
	
	
	public LargeHashMap(String libraryObjectFile){
		this.libraryObjectFile = libraryObjectFile;
		this.table = new HashMap();
		this.keys = new TreeSet();
		this.beginTable = new  HashMap();
		this.endTable = new HashMap();
		try{
			this.RAF = new RandomAccessFile(this.libraryObjectFile, "rw");
			index1 = 0;//RAF.getFilePointer();
			//RAF.writeLong(1L);  
			index2 = 8;//RAF.getFilePointer();
			//RAF.writeLong(1L);
			index3 = 16;//RAF.getFilePointer();
			//RAF.writeLong(1L);  
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		//this.loadLibraryFromFile(this.libraryObjectFile);
	}
	
	
	public void loadLibraryFromFile(String file){
		this.libraryObjectFile = file;
		try{
			if(this.RAF != null){
				this.RAF.close();
			}
			this.RAF = new RandomAccessFile(file, "r");
			index1 = RAF.readLong();
			index2 = RAF.readLong();
			index3 = RAF.readLong();
			long current = RAF.getFilePointer();
			if(DEBUG){
				System.out.println("index is: " + index1 + "\t" + index2 + "\t" + index3);
			}
			RAF.seek(index1);
			byte[] input = new byte[(int)(index2-index1)];
			if(DEBUG){
				System.out.println("current is : " + RAF.getFilePointer());
			}
			int readed = RAF.read(input);
			if(DEBUG){
				System.out.println("read in bytes: " + readed);
			}
			ByteArrayInputStream in = new ByteArrayInputStream(input);
			ObjectInputStream in2 = new ObjectInputStream(in);
			this.beginTable = (Map)in2.readObject();
			in.close();
			in2.close();
			RAF.seek(index2);
			input = new byte[(int)(index3-index2)];
			if(DEBUG){
				System.out.println("current is : " + RAF.getFilePointer());
			}
			readed = RAF.read(input);
			if(DEBUG){
				System.out.println("read in bytes: " + readed);
			}
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
			if(this.RAF == null){
				this.RAF = new RandomAccessFile(this.libraryObjectFile, "rw");
			}
			//first we store file position to beginTable and endTable
			ByteArrayOutputStream out1;
			ObjectOutputStream out2;
			int counter=0;
			for(Iterator it = specTable.keySet().iterator(); it.hasNext();){
				Object key = it.next();
				//beginTable.put(key, RAF.getFilePointer());
				//System.out.println("begin is: " + RAF.getFilePointer());
				Object values = specTable.get(key);
				this.put(key, values);
				//out1 = new ByteArrayOutputStream();
				//out2 = new ObjectOutputStream(out1);
				//out2.writeObject(values);
				//out2.flush();
				//RAF.write(out1.toByteArray());
				//endTable.put(key, RAF.getFilePointer());
				//System.out.println("end is : " + RAF.getFilePointer());
			}
			//storing the index table
			finalize();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public void finalize(){
		try{
			ByteArrayOutputStream out1;
			ObjectOutputStream out2;
			long current = RAF.getFilePointer();
			RAF.seek(index1);
			RAF.writeLong(current);
			RAF.seek(current);
			if(DEBUG){
				System.out.println("begin writing " + "@" + index1 + ": " + current);
			}
			out1 = new ByteArrayOutputStream();
			out2 = new ObjectOutputStream(out1);
			out2.writeObject(this.beginTable);
			out2.flush();
			RAF.write(out1.toByteArray());
			if(DEBUG){
				System.out.println("end writing " + RAF.getFilePointer() + " supposely written: " + out1.toByteArray().length);
			}
			out1.close();
			out2.close();
		
			current = RAF.getFilePointer();
			RAF.seek(index2);
			RAF.writeLong(current);
			RAF.seek(current);
			if(DEBUG){
				System.out.println("begin writing " + "@" + index2 + ":" + current);
			}
			out1 = new ByteArrayOutputStream();
			out2 = new ObjectOutputStream(out1);
			out2.writeObject(this.endTable);
			out2.flush();	
			RAF.write(out1.toByteArray());
			if(DEBUG){
				System.out.println("end writing " + RAF.getFilePointer() + " supposely written: " + out1.toByteArray().length);
			}
			current = RAF.getFilePointer();
			if(DEBUG){
				System.out.println("begin writing " + "@" + index3 + ":" + current);
			}
			RAF.seek(index3);
			RAF.writeLong(current);
			RAF.seek(current);
			out1.close();
			out2.close();
			out2.writeObject(this.beginTable);
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
	
	public void put(Object key, Object values){
		try{
			ByteArrayOutputStream out1;
			ObjectOutputStream out2;
			if(this.lastEntryIndex==0){
				RAF.writeLong(1L);
				RAF.writeLong(1L);
				RAF.writeLong(1L);
			}
			beginTable.put(key, RAF.getFilePointer());
			//System.out.println("begin is: " + RAF.getFilePointer());
			out1 = new ByteArrayOutputStream();
			out2 = new ObjectOutputStream(out1);
			out2.writeObject(values);
			out2.flush();
			RAF.write(out1.toByteArray());
			endTable.put(key, RAF.getFilePointer());
			this.lastEntryIndex = RAF.getFilePointer();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
	
	public Set getKeys(){
		return this.beginTable.keySet();
	}

	public String getLibraryObjectFile() {
		return libraryObjectFile;
	}


	public void setLibraryObjectFile(String libraryObjectFile) {
		this.libraryObjectFile = libraryObjectFile;
	}

}
