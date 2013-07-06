package MSPLIT.SpectrumIO;
import MSPLIT.Spectrum;

public class MGFParser extends SpectrumParser{

	public MGFParser(String filename) {
		super(filename);
	}

	@Override
	protected Spectrum getNextSpectrum() {
		return SpectrumReader.readSpectrumFromMGF(this.getBuff());
	}
	
}
