package org.Spectrums;

/**
 * Contain Utiliy function for array
 * @author jian wang
 *
 */
public class ArrayUtils {
	
	public static void addArray(double[] sum, double[] toAdd){
		for(int i = 0; i < sum.length; i++){
			sum[i] = sum[i] + toAdd[i];
		}
	}

	
	public static void divideArray(double[] arry, double factor){
		for(int i = 0; i < arry.length; i++){
			arry[i] /= factor;
		}
	}
	
	public static void divideArray(double[][] arry, double factor){
		for(int i = 0; i < arry.length; i++){
			for(int j = 0; j < arry[i].length; j++){
				arry[i][j] /= factor;
			}
		}
	}
	
	public static void normalizeArray(double[] arry){
		divideArray(arry, sum(arry));
	}
	
	public static double sum(double[] arry){
		double sum = 0.0;
		for(int i = 0; i < arry.length; i++){
			sum += arry[i];
		}
		return sum;
	}
	
	public static double sum(double[][] arry){
		double total = 0.0;
		for(int i = 0; i < arry.length; i++){
			total += sum(arry[i]);
		}
		return total;
	}
	
	public static double[] shift(double[] arry, double shift){
		for(int i = 0; i < arry.length; i++){
			arry[i] += shift;
		}
		return arry;
	}
	
	//added a very small number to avoid div-by-zero errors
	public static void addPseudoCounts(double[][] counts){
		for(int i = 0; i < counts.length; i++){
			for(int j = 0; j < counts[i].length; j++){
				counts[i][j] += counts[i][j] == 0 ? 0.1 : 0; 
			}
		}
	}
	
	public static double[] mergeArray(double[] array1, double[] array2){
		double[] combine = new double[array1.length + array2.length];
		for(int i = 0, size = array1.length; i < size; i++){
			combine[i] = array1[i];
		}
		for(int i = 0, offset= array1.length, size = array2.length; 
			i < size; i++){
			combine[i+offset] = array2[i];
		}
		return combine;
	}
	
	/**
	 * Given an element and a array of elements which represent the
	 * boundaries of intervals, return the interval that this element
	 * fall into, noted by default interval is left inclusive
	 * @param element
	 * @param intervals
	 * @return
	 */
	public static int getIntervalIndex(int element, int[] intervals){
		for(int i = 1, size = intervals.length; i<size; i++){
			if(element < intervals[i] && element >= intervals[i-1]){
				return i-1;
			}
		}
		throw new IllegalArgumentException("element not within Interval: " + element);
	}
	
	public static int getIntervalIndex(double element, double[] intervals){
		for(int i = 1, size = intervals.length; i<size; i++){
			if(element < intervals[i] && element >= intervals[i-1]){
				return i-1;
			}
		}
		throw new IllegalArgumentException("element not within Interval: " + element);
	}
	
	public static void printTable(double[][] arry){
		for(int i = 0; i < arry.length; i++){
			for(int j = 0; j < arry[i].length; j++){
				System.out.print(arry[i][j] + "\t");
			}
			System.out.println();
		}
	}
}
