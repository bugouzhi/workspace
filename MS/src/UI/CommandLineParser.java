package UI;

import java.util.HashMap;
import java.util.Map;

import Utils.StringUtils;

/**
 * This class process unix-like command line arguments 
 * @author Jian
 *
 */
public class CommandLineParser {
	String[] args;
	Map<String, String> argValuePairs;
	public CommandLineParser(String[] args){
		this.args = args;
		this.argValuePairs = new HashMap();
		if(this.args.length % 2 > 0){
			//System.err.println("Command line arguments and value pairs do not match...");
			//return;
		}else{
			for(int i = 0; i < args.length; i=i+2){
				this.argValuePairs.put(args[i], args[i+1]);
			}
		}
	}
	
	public double getDouble(String option){
		if(argValuePairs.containsKey(option)){
			return Double.parseDouble(argValuePairs.get(option));
		}else{
			throw new IllegalArgumentException("Option: " + option + " is not specified");
		}
	}
	
	
	public int getInteger(String option){
		if(argValuePairs.containsKey(option)){
			return Integer.parseInt(argValuePairs.get(option));
		}else{
			throw new IllegalArgumentException("Option: " + option + " is not specified");
		}
	}
	
	public String getString(String option){
		if(argValuePairs.containsKey(option)){
			return argValuePairs.get(option);
		}else{
			throw new IllegalArgumentException("Option: " + option + " is not specified");
		}
	}
	
	public double getDouble(int argumentInd){
		if(argumentInd < this.args.length){
			return Double.parseDouble(this.args[argumentInd]);
		}else{
			throw new IllegalArgumentException("There is no argument at" + argumentInd + "position");
		}
	}
	
	
	public int getInteger(int argumentInd){
		if(argumentInd < this.args.length){
			return Integer.parseInt((this.args[argumentInd]));
		}else{
			throw new IllegalArgumentException("There is no argument at" + argumentInd + "position");
		}
	}
	
	public String getString(int argumentInd){
		if(argumentInd < this.args.length){
			return (this.args[argumentInd]);
		}else{
			throw new IllegalArgumentException("There is no argument at" + argumentInd + "position");
		}
	}
	
	public int[] getIntegers(int argumentInd){
		if(argumentInd < this.args.length){
			String[] tokens = args[argumentInd].split(",");
			int[] values = StringUtils.toIntArray(tokens);
			return values;
		}else{
			throw new IllegalArgumentException("There is no argument at" + argumentInd + "position");
		}
	}
	
	public boolean isAvailable(int ind){
		return this.args.length > ind;
	}
	public double[] getDoubles(int argumentInd){
		if(argumentInd < this.args.length){
			String[] tokens = args[argumentInd].split(",");
			double[] values = StringUtils.toDoubleArray(tokens);
			return values;
		}else{
			throw new IllegalArgumentException("There is no argument at" + argumentInd + "position");
		}
	}
	
}
