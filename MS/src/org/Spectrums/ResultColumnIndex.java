package org.Spectrums;
/**
 * A class search result column index used to parse result file
 * @author Jian
 *
 */
public class ResultColumnIndex {
	
		private int specInd;
		private int pepInd;
		private int protInd;
		private int fileInd;
		private int scoreInd;
		private int chargeInd;
		private int sortOrder;
		public static ResultColumnIndex MSGFDB_INDEX = new ResultColumnIndex(0, 1, 7, 8, 6, 11, 1);
		//public static ResultColumnIndex SWATHMSPLIT_INDEX = new ResultColumnIndex(0, 1, 7, 8, 6, 12, 1);
		
		public ResultColumnIndex(int fileInd, int specInd, int pepInd, int protInd, int chargeInd, int scoreInd, int sortOrder){
			this.fileInd = fileInd;
			this.specInd = specInd;
			this.pepInd = pepInd;
			this.protInd = protInd;
			this.chargeInd = chargeInd;
			this.scoreInd  = scoreInd;
			this.sortOrder = sortOrder;
		}
		
		
		public int getSortOrder() {
			return sortOrder;
		}

		public void setSortOrder(int sortOrder) {
			this.sortOrder = sortOrder;
		}
		
		public int getSpecInd() {
			return specInd;
		}
		public void setSpecInd(int specInd) {
			this.specInd = specInd;
		}
		public int getPepInd() {
			return pepInd;
		}
		public void setPepInd(int pepInd) {
			this.pepInd = pepInd;
		}
		public int getProtInd() {
			return protInd;
		}
		public void setProtInd(int protInd) {
			this.protInd = protInd;
		}
		public int getFileInd() {
			return fileInd;
		}
		public void setFileInd(int fileInd) {
			this.fileInd = fileInd;
		}
		public int getScoreInd() {
			return scoreInd;
		}
		public void setScoreInd(int scoreInd) {
			this.scoreInd = scoreInd;
		}
		public int getChargeInd() {
			return chargeInd;
		}
		public void setChargeInd(int chargeInd) {
			this.chargeInd = chargeInd;
		}
		

}
