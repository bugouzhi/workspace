package Spectrum;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

/**
 * Mass information for amino acid and various type of ions
 * @author jian wang
 *
 */
public class Mass {
	public static double PROTON_MASS = 1.007276;
	public static double DSPLINKER_MASS = 138.0680;//173.9809;
	public static double DEUTERIUM = 2.0136;
	public static double C13 = 13.00335;
	public static double C12 = 12.0000;
	public static double WATER = 18.0106;
	public static double NH3 = 17.0265;
	public int AA_ALPHABET_SIZE = 30;
	public static HashMap<String, Double> modMap = initialize();
	public static double[] aaMap = initializeAA();
	public static String[] standardIonsType = {"b", "b-H20", "b-NH3", "b(iso)", "b-H20-H20", "b-H20-NH3", "a", "a-H20", "a-NH3",
	"y", "y-H20", "y-NH3", "y(iso)", "y-H20-H20", "y-H20-NH3", "Noise"};
	public static String[] standardPrefixes = {"b", "b-H20", "b-NH3", "b(iso)", "b-H20-H20", "b-H20-NH3"};
	public static String[] standardSuffixes = {	"y", "y-H20", "y-NH3", "y(iso)", "y-H20-H20", "y-H20-NH3"};
	public static String[] standardPrefixesX = {"b",  "b(iso)", "b(X)", "b(Xiso)", "b-H20", "b-NH3","b-H20-H20", "b-H20-NH3"}; //with isotop-coded xlinker
	public static String[] standardSuffixesX = {"y", "y(iso)", "y(X)", "y(Xiso)", "y-H20", "y-NH3","y-H20-H20", "y-H20-NH3"};

	public static double maxAAMass = 200;	
	
	/**
	 * Load masses information for amino acid and Ions from resources files
	 * @return
	 */
	public static HashMap initialize(){
		String file = "/Resources/IonsMod.txt";
		modMap = new HashMap();
		try{
			InputStream in = Mass.class.getResourceAsStream(file);
			if(in == null){
				System.out.println("cannot find resources");
				return null;
			}
			BufferedReader bf = new BufferedReader(new InputStreamReader(in));
			String line = bf.readLine();
			String[] tokens;
			while(line != null){
				tokens = line.split("\t");
				modMap.put(tokens[0], new Double(tokens[1]));
				line = bf.readLine();
			}
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
		return modMap;
	}
	
	/**
	 * Initialize amino acid masses table
	 * @return
	 */
	public static double[] initializeAA(){
		String file = "/Resources/AAmass.txt";
		aaMap = new double[30];
		try{
			InputStream in = Mass.class.getResourceAsStream(file);
			if(in == null){
				System.out.println("cannot find resources");
				return null;
			}
			BufferedReader bf = new BufferedReader(new InputStreamReader(in));
			String line = bf.readLine();
			String[] tokens;
			while(line != null){
				tokens = line.split("\t");
				aaMap[(int)(tokens[0].charAt(0)-'A')] = Double.parseDouble(tokens[1]); 
				line = bf.readLine();
			}
		}catch(IOException ioe){
			System.out.println(ioe.getMessage());
			ioe.printStackTrace();
		}
//		for(int i = 0; i < aaMap.length; i++){
//			if(aaMap[i] == 0){
//				aaMap[i] = 100000;              //set unrecognize aa alphabet to huge masses, so we never consider such peptide
//			}
//		}
		return aaMap;
	}
	
	/**
	 * Give ion type return the mass offset
	 * @param ion
	 * @return
	 */
	public static double getIonMod(String ion){
		if(modMap.containsKey(ion)){
			return modMap.get(ion).doubleValue();
		}else{
			System.err.println("warning: ion type do not exits, returning zero mass mod");
			return 0.0;
		}
	}
	
	/**
	 * 
	 * @param aa
	 * @return
	 */
	public static double getAAMass(char aa){
		//System.out.println("aa: " + aa);
		if(aa == '*'){
			return 100000;
		}
		if(aa == 'U'){
			return 100000;
		}
		if((int)aa-'A' < 0){
			System.out.println("warning: " + aa + " might not be a valide aa");
		}
		return aaMap[(int)aa-'A'];
	}
	
	/**
	 * check if mass is a valide aa mass
	 * @param mass
	 * @return
	 */
	public static boolean isAAMass(double mass){
		for(int i = 0; i < Mass.aaMap.length; i++){
			if(Math.abs(mass - aaMap[i]) < 1.0){
				return true;
			}
		}
		return false;
	}
	
	public static void main(String[] args){
		System.out.println("hihi");
	}
	
	/**
	 * 
	 * @param mass
	 * @param decimal
	 * @return
	 */
	public static double roundMass(double mass, int decimal){
		System.out.println("started: " + mass);
		long round = Math.round((mass*10*decimal));
		System.out.println("long is: " + round);
		return ((double)round)/(10*decimal);
	}
	
	/**
	 * Compute mass difference in either Da mode or ppm mode
	 * @param mass1
	 * @param mass2
	 * @param mode
	 * @return
	 */
	public static double massDiff(double mass1, double mass2, int mode){
		if(mode == 1){
			return Math.abs(mass1 - mass2);
		}else{
			return Math.abs((mass1 - mass2)*1000000 / mass2); 
		}

	}
	
}
