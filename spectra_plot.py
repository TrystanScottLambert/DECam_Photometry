""" Performing synthetic photometry on the spectrum of the quasar."""

import numpy as np
import astropy.units as u
import pylab as plt
from synphot import Observation, SpectralElement, SourceSpectrum
from synphot.models import Empirical1D



if __name__ == '__main__':
    SPECTRA_FILE = '../QSO_Spectra/J2348_dered.txt'
    BANDPASS_FILE = '../QSO_Spectra/NB964_DECam_F29.txt'
    band_pass = SpectralElement.from_file(BANDPASS_FILE, wave_unit='nm')
    wave, flux, err = np.loadtxt(SPECTRA_FILE, unpack=True)
    wave = wave * u.angstrom
    flux = flux * (1e-17 * u.erg)/u.s/(u.cm**2)/u.angstrom
    spectrum = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)

    observation = Observation(spectrum, band_pass)
    print(observation.effstim(u.ABmag))
