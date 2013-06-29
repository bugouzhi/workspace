package mixdb;
/**
 * Test suffix array
 * @author Jian Wang
 *
 */
import msdbsearch.*;
import sequences.*;

public class SuffArrayTest {
	public static void testGetSuffArray(String dbPath){
		dbPath = "../mixture_linked/database/Human_allproteins.fasta";
		CompactFastaSequence sequence = new CompactFastaSequence(dbPath);
		CompactSuffixArray suffarry = new CompactSuffixArray(sequence);
		sequence = suffarry.getSequence();
		System.out.println("number of proteins: " + sequence.getNumProteins());
		System.out.println("suffixes siez: " + suffarry.getSize());
		for(long i = 0; i <= sequence.getSize(); i++){
			//System.out.print(sequence.getCharAt(i));
		}
	}
	
	public static void testFastaSeq(String dbPath){
		dbPath = "../mixture_linked/database/Human_allproteins.fasta";
		FastaSequence sequence = new FastaSequence(dbPath);
		//CompactFastaSequence sequence = new CompactFastaSequence(dbPath);
		System.out.println("size: " + sequence.getSize());
		long size = sequence.getSize();
		for(int i = 0; i < 30; i++){
		for(long j = 0; j < size; j++){
			if(i % 1000 == 0){
				//System.out.println();
				//System.out.println(i);
			}
			//System.out.print(sequence.getCharAt(i));
		}
		System.out.println(i);
		}
	}
	
	public static void main(String[] args){
		//testGetSuffArray("");
		testFastaSeq("");
	}
}
