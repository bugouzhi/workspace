package IO;
import java.util.Iterator;
import Spectrum.Spectrum;
/**
 * Generaric interface of spectrum reader
 * @author Jian
 *
 */
public interface SpectrumReader extends Iterator<Spectrum>{
	public Spectrum readSpectrumByIndex(int index);
}
