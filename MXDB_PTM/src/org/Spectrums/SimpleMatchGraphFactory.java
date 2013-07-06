package org.Spectrums;

import java.util.ArrayList;
import java.util.List;

public class SimpleMatchGraphFactory {
	private static int currentIndex = 0;
	private static List<SimpleMatchingGraph> graphList = new ArrayList<SimpleMatchingGraph>();
	private static SimpleMatchingGraph empty = new SimpleMatchingGraph();
	
	public static SimpleMatchingGraph createSimpleMatchGraph(){
		if(currentIndex < graphList.size()){
			SimpleMatchingGraph g = graphList.get(currentIndex);
			g.clearGraph();
			return g;
		}else{
			SimpleMatchingGraph g = new SimpleMatchingGraph();
			graphList.add(g);
			return g;
		}
		
	}
	
	public static void resetFactory(){
		currentIndex = 0;
	}
	
	public static void clearFactory(){
		graphList.clear();
	}
	
	public static SimpleMatchingGraph getEmptyGraph(){
		return empty;
	}
}
