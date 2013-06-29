package mixdb;

import org.Spectrums.PeptideLite;

/**
 * An array based implementation of the theoretical spectrum
 * allow faster process
 * @author Jian Wang
 *
 */
public class ArrayTheoreticalSpectrum extends ArraySpectrum{
	private PeptideLite peplite;

	public PeptideLite getPeplite() {
		return peplite;
	}

	public void setPeplite(PeptideLite peplite) {
		this.peplite = peplite;
	}
	
}
