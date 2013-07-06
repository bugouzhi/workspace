package MSPLIT.SpectrumIO;

import MSPLIT.Spectrum;

public class MS2Parser extends SpectrumParser{

	public MS2Parser(String filename) {
		super(filename);
	}

	@Override
	protected Spectrum getNextSpectrum() {
		return SpectrumReader.readSpectrumFromMS2(this.getBuff());
	}

}
