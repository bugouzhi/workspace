package org.Spectrums;
import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.NativeLibrary;
import com.sun.jna.Pointer;

public class WiffSDKAccess {
	public static String wiffSDKPath = "C:\\Documents and Settings\\Jian Wang\\workspace\\WiffSDK\\WiffReaderSDK\\Bin";
	
	public interface WiffReader extends Library{
		
	}
	
	public static void main(String[] args){
		System.setProperty("jna.library.path", wiffSDKPath);
		WiffReader wiff = (WiffReader)Native.loadLibrary("Clearcore2.Data", WiffReader.class);
		NativeLibrary wiffNative = NativeLibrary.getInstance("Clearcore2.Data");
		//Pointer license = wiffNative.getGlobalVariableAddress("Licensing");
		//String[] keys = license.getStringArray(0);
		//for(int i = 0; i < keys.length; i++){
		//	System.out.println(keys[i]);
		//}
	}
}
