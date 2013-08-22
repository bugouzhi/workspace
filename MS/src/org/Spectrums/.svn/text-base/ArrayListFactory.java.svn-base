package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

public class ArrayListFactory {
	private static int currentIndex = 0;
	private static List<ArrayList> lists = new ArrayList<ArrayList>();
	
	public static List createList(){
		if(currentIndex < lists.size()){
			ArrayList l = lists.get(currentIndex);
			l.clear();
			return l;
		}else{
			ArrayList l = new ArrayList();
			lists.add(l);
			return l;
		}
		
	}
	
	public static void resetFactory(){
		currentIndex = 0;
	}
	
	public static void clearFactory(){
		lists.clear();
	}
}
