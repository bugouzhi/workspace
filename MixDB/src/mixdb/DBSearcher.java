package mixdb;

import java.util.Collection;
import java.util.TreeSet;

import org.Spectrums.Spectrum;

/**
 * A general interface for searching a spectrum
 * against a database
 * @author Jian Wang
 *
 */
public interface DBSearcher {
	public TreeSet topNCandidates(Spectrum s, Collection db, int topN);
}
