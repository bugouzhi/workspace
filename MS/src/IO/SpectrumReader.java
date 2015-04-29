package IO;
import java.util.Iterator;

import org.Spectrums.Spectrum;

public interface SpectrumReader extends Iterator<Spectrum>{
	public Spectrum readSpectrumByIndex(int index);
}
