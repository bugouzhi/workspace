package UI;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * This class run some command line using the shell
 * @author Jian
 *
 */
public class CommandLineRunner {
	public static void runCommand(String cmd){
		try{
			System.out.println(cmd);
			String[] commands = new String[]{cmd};
			Process p1 = Runtime.getRuntime().exec(commands);
			gobbleStream(p1.getInputStream());
			gobbleStream(p1.getErrorStream());
			p1.waitFor();
		}catch(Exception e){
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
	}
	
	//need to extract output from svm_classify otherwise process will not terminate
	private static void gobbleStream(InputStream is){
		try{
			BufferedReader buff =new BufferedReader(new InputStreamReader(is));
			String line = buff.readLine();
			while(line != null){
				System.out.println(line);
				line = buff.readLine();
			}
		}catch(IOException ioe){
			System.err.println(ioe.getMessage());
			ioe.printStackTrace();
		}
	}
}
