package MSPLIT.SpectrumIO;

import MSPLIT.Spectrum;

public class MSPParser extends SpectrumParser{

	public MSPParser(String filename) {
		super(filename);
	}

	@Override
	protected Spectrum getNextSpectrum() {
		return SpectrumReader.readSpectrumFromMSP(this.getBuff());
	}
	
}
