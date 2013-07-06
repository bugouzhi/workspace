package MSPLIT.SpectrumIO;

import MSPLIT.Spectrum;

public class SplibParser extends SpectrumParser{

	public SplibParser(String filename) {
		super(filename);
	}

	@Override
	protected Spectrum getNextSpectrum() {
		return SpectrumReader.readSpectrumFromSplib(this.getBuff());
	}
	
}
