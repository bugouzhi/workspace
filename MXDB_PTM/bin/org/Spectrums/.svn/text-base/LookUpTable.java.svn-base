package org.Spectrums;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;

/**
 * General N-dimensional table, backed by an one dimensional array
 * access to the table through index manipulation
 * @author Jian Wang
 *
 */
public class LookUpTable implements Serializable{
	/**
	 * for Serializable interface
	 */
	private static final long serialVersionUID = 5699911929597533852L;
	private double[] table;
	private int[] dim;
	public LookUpTable(int[] dimension){
		this.dim = dimension;
		int size = 1;
		for(int i = 0; i < dimension.length; i++){
			size *= dimension[i];
		}
		this.table = new double[size];
	}
	
	public double get(int[] index){
		int realIndex = getIndex(index);
		//System.out.println("index are : " + index[0] + "," + index[1] + "," + index[2]);
		if(realIndex >= this.table.length){
			System.err.print("Index exceed table dimension: ");
			for(int i = 0; i < index.length; i++){
				System.err.print(index[i] + "\t");
			}
			System.err.println();
			throw new IllegalArgumentException();
		}
		return table[realIndex];
	}
	
	public void put(int[] index, double value){
		int realIndex = getIndex(index);
		table[realIndex] = value;
	}
	
	private int getIndex(int[] index){
		if(index.length != dim.length){
			throw new IllegalArgumentException("Dimension of index and table must match");
		}
		int realIndex = 0;
		int scale = 1;
		for(int i = index.length-1; i >= 0; i--){
			realIndex +=  index[i]*scale;
			scale *= this.dim[i];
			//System.out.println("real: " + realIndex + "\t" + "scale: " + scale);
		}
		return realIndex;
	}
	
	/**
	 * increment the count in the particular entry by one
	 * @param lp
	 * @param table
	 */
	public void incrementIonCount(int[] index){
		double currCount = this.get(index);
		//System.out.println("before increment: " + currCount);
		this.put(index, currCount+1);
//		System.out.println("after increment: " + this.get(index));
	}
	
	
	public void incrementIonCount(int[] index, double increment){
		double currCount = this.get(index);
		//System.out.println("before increment: " + currCount);
		this.put(index, currCount+increment);
	}
	
	public static void mapOperator(LookUpTable table1, LookUpTable table2, TableElementOperator op){
		if(table1.table.length != table2.table.length){
			throw new IllegalArgumentException("table for map operator must be or same dimension");
		}
		
		for(int i = 0; i < table1.table.length; i++){
			table1.table[i] = op.eval(table1.table[i], table2.table[i]);
		}
		
	}
	
	public static void mapOperator(LookUpTable table1, double val, TableElementOperator op){
		for(int i = 0; i < table1.table.length; i++){
			table1.table[i] = op.eval(table1.table[i], val);
		}
		
	}
	
	public void writeLibToFile(String outfile){
		try{
			BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(outfile));
			ObjectOutputStream oo = new ObjectOutputStream(bo);
		    oo.writeObject(this);
		    oo.flush();
		    oo.close();
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			
		}
	}
}
