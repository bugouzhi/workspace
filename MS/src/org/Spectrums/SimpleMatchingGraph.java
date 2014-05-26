package org.Spectrums;
/**
 * Very basic implementation of a bi-partite graph object
 * using adacency list. We try to keep the implementation
 * minimal as to optimize the performance for mathicng graph
 * @author jian wang
 *
 */
import java.util.Collection;
import java.util.LinkedList;
import java.util.Map;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;
import java.util.Iterator;
import org.jgrapht.EdgeFactory;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;

public class SimpleMatchingGraph implements UndirectedGraph{
	public static int Observed = 1;
	public static int Theoretical = 2;
	private Map<Object, List> adjList1; //by convention group1 is actual peaks
	private Map<Object, List> adjList2; //on the other hand group2 is theoretical peaks
	private List<Object> list1;
	
	public SimpleMatchingGraph(){
		this.adjList1 = new HashMap(200);
		this.adjList2 = new HashMap(200);
		
	}
	
	public SimpleMatchingGraph(SimpleMatchingGraph g1, SimpleMatchingGraph g2){
		this();
		copyTable(g1.adjList2, this.adjList2);
		copyTable(g2.adjList2, this.adjList2);
		copyTable(g1.adjList1, this.adjList1);
		copyTable(g2.adjList1, this.adjList1);
	}
	
	//copy if not exists, append if they do
	private void copyTable(Map source, Map target){
		Set vertexSet = source.keySet();
		Iterator it = vertexSet.iterator();
		Object key;
		List l;
		while(it.hasNext()){
			key = it.next();
			if(target.containsKey(key)){
				l = (List)target.get(key);
			}else{
				l = new ArrayList();
			}
			l.addAll((List)source.get(key));
			target.put(key, l);
		}
	}
	
	public boolean addVertex(Object o, int group){
		if(group == 1){
			if(!adjList1.containsKey(o)){
				adjList1.put(o, new ArrayList());
				//adjList1.put(o, ArrayListFactory.createList());
				
			}
		}else{
			if(!adjList2.containsKey(o)){
				adjList2.put(o, new ArrayList());
				//adjList2.put(o, ArrayListFactory.createList());
				
			}
		}
		return true;
	}
	
	public boolean addVertex(Object o){
		if(!adjList1.containsKey(o)){
			adjList1.put(o, new ArrayList());
		}
		return true;
	}
	
	/**
	 * first vertex in group A, second vertex in group B
	 */
	public Object addEdge(Object v1, Object v2){
		if(adjList1.containsKey(v1) && adjList2.containsKey(v2)){
			adjList1.get(v1).add(v2);
			adjList2.get(v2).add(v1);
		}else{
			
			//System.out.println("vertex not exists in graph");
		}
		return new DefaultEdge();
	}
	
	public List getNeighbors(Object v){
		if(this.adjList1.containsKey(v)){
			return adjList1.get(v);
		}else{
			return adjList2.get(v);
		}	
	}
	
	public void removeAllNeighbor(Object v){
		List<Peak> neigh;
		if(this.adjList1.containsKey(v)){
			neigh = adjList1.get(v);
		}else{
			neigh = adjList2.get(v);
		}
		for(int i = 0; i < neigh.size(); i++){
			Peak currentNeigh = neigh.get(i);
			List<Peak> neigh2 = this.getNeighbors(currentNeigh);
			neigh2.remove(v);
		}
		neigh.clear();
	}
	
	public Collection getVerticeWithEdges(int group, int minEdges){
		Set vertexSet = new HashSet();
		Iterator it;
		Map<Object, List> adjList;
		if(group == SimpleMatchingGraph.Observed){
			adjList = this.adjList1;
		}else{
			adjList = this.adjList2;
		}
		for(it = adjList.keySet().iterator(); it.hasNext();){
			Object key = it.next();
			if(adjList.containsKey(key) && adjList.get(key).size() >= minEdges){
				vertexSet.add(key);
			}
		}
		return vertexSet;
	}
	
	/**
	 * create a bi-partite graph with each node
	 * has at most one edge, use closest mass to determine
	 * which peaks got matched
	 */
	public void toBiPartiteByMassError(){
		for(Iterator<LabelledPeak> it = this.vertexSet(SimpleMatchingGraph.Observed).iterator(); it.hasNext();){
			Peak p = it.next();
			List neighs = this.getNeighbors(p);
			LabelledPeak closestP = null;
			double smallestErr = 10000;
			for(int i = 0; i < neighs.size(); i++){
				LabelledPeak currPeak = (LabelledPeak)neighs.get(i);
				double currDiff = Math.abs(currPeak.getMass()-p.getMass());
				closestP = currDiff < smallestErr ? currPeak : closestP;
				smallestErr = currDiff < smallestErr ? currDiff : smallestErr;
			}
			if(closestP != null){
				neighs.clear();
				neighs.add(closestP);
			}
		}
	}
	
	
	public void clearGraph(){
		this.adjList1.clear();
		this.adjList2.clear();
	}
	@Override
	public boolean addEdge(Object sourceVertex, Object targetVertex, Object e) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean containsEdge(Object sourceVertex, Object targetVertex) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean containsEdge(Object e) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean containsVertex(Object v) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set edgeSet() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Set edgesOf(Object vertex) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public Set getAllEdges(Object sourceVertex, Object targetVertex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getEdge(Object sourceVertex, Object targetVertex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public EdgeFactory getEdgeFactory() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getEdgeSource(Object e) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getEdgeTarget(Object e) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getEdgeWeight(Object e) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public boolean removeAllEdges(Collection edges) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set removeAllEdges(Object sourceVertex, Object targetVertex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean removeAllVertices(Collection vertices) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Object removeEdge(Object sourceVertex, Object targetVertex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean removeEdge(Object e) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean removeVertex(Object v) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set vertexSet() {
		Set s = new HashSet();
		s.addAll(adjList1.keySet());
		s.addAll(adjList2.keySet());
		return s;
	}
	
	public Set vertexSet(int group){
		if(group == 1){
			return this.adjList1.keySet();
		}else{
			return this.adjList2.keySet();
		}
	}

	@Override
	public int degreeOf(Object vertex) {
		List neighbors; 
		if(this.adjList1.containsKey(vertex)){
			neighbors = this.adjList1.get(vertex);
		}else{
			neighbors = this.adjList2.get(vertex);
		}
		if(neighbors == null){
			return 0;
		}else{
			return neighbors.size();
		}
	}
	
	public Set getNeighborSet(Object vertex){
		List l = null;
		Set neighbors = new HashSet();
		if(this.adjList1.containsKey(vertex)){
			l = this.adjList1.get(vertex);
		}else{
			l = this.adjList2.get(vertex);
		}
		neighbors.addAll(l);
		return neighbors;
	}
	
}
