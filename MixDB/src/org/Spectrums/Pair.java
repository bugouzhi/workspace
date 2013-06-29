package org.Spectrums;
/**
 * Generaic class for storing a pair of object
 * @author Jian Wang
 *
 */
public class Pair {
	private Object first;
	private Object second;
	public Pair(Object o1, Object o2){
		this.first = o1;
		this.second = o2;
	}
	public Object getFirst() {
		return first;
	}
	public void setFirst(Object first) {
		this.first = first;
	}
	public Object getSecond() {
		return second;
	}
	public void setSecond(Object second) {
		this.second = second;
	}
	
	public boolean equals(Object o){
		Pair p = (Pair)o;
		return (p.first.equals(this.getFirst()) && p.second.equals(this.getSecond()))
			|| (p.first.equals(this.getSecond()) && p.second.equals(this.getFirst()));
	}
	
	//need to reimplement this
	public int hashCode(){
		return first.hashCode() + second.hashCode();
	}
	
	
}
