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
	
	//normalize by sum of all values
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
	
	public static String getString(double[] arry){
		return getString(arry, "%1$.2f"+"\t");
	}
	public static String getString(double[] arry, String format){
		StringBuffer str = new StringBuffer();
		for(int i = 0; i < arry.length; i++){
			//str.append(String.format("%1$.2f\t", arry[i]));
			str.append(String.format(format+"\t", arry[i]));
		}
		return str.toString();
	}
	
	public static double dotProd(double[] arry1, double[] arry2){
		double sum = 0;
		for(int i = 0; i < arry1.length; i++){
			sum+=arry1[i]*arry2[i];
		}
		return sum;
	}
	
	public static double average(double[] arry1){
		double sum = 0;
		for(int i = 0; i < arry1.length; i++){
			sum+=arry1[i];
		}
		//return 0;
		return sum/((double)arry1.length);
	}
	
	public static double[] offSet(double[] arry1, double offset, double[] arry2){
		for(int i = 0; i < arry1.length; i++){
			arry2[i] = arry1[i] + offset;
		}
		return arry2;
	}
	
	//normalize by euclidean norm
	public static void normalize(double[] arry){
		double mag = 0;
		for(int j = 0; j < arry.length; j++){
			mag += arry[j]*arry[j];
		}
		mag = Math.pow(mag, 0.5);
		mag += 0.000000000001;//avoid div by zero
		for(int j = 0; j < arry.length; j++){
			arry[j] /= mag;
		}	
	}
	
	public static void sqrt(double[] arry){
		for(int j = 0; j < arry.length; j++){
			if(arry[j] > 0){
				arry[j] = Math.pow(arry[j], 0.5);
			}else{
				arry[j] = 0;
			}
		}	
	}
	
	public static double[] getArray(Spectrum s){
		//convert lib spect into array
		double[] sarry = new double[s.getPeak().size()];
		for(int i = 0; i < s.getPeak().size(); i++){
				sarry[i] = s.getPeak().get(i).getIntensity();
		}
		return sarry;
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
