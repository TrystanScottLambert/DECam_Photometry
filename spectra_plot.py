""" Performing synthetic photometry on the spectrum of the quasar."""

import numpy as np
import astropy.units as u
from synphot import Observation, SpectralElement, SourceSpectrum
from synphot.models import Empirical1D



if __name__ == '__main__':
    SPECTRA_FILE = '../QSO_Spectra/J2348_dered.txt'
    NARROW_BANDPASS_FILE = '../QSO_Spectra/NB964_DECam_F29.txt'
    I_BANDPASS_FILE = '../QSO_Spectra/decam_i_bandpass.txt'
    Z_BANDPASS_FILE = '../QSO_Spectra/decam_z_bandpass.txt'


    #Define the filter objects
    narrow_band_pass = SpectralElement.from_file(NARROW_BANDPASS_FILE, wave_unit='nm')
    i_band_pass = SpectralElement.from_file(I_BANDPASS_FILE, wave_unit='nm')
    z_band_pass = SpectralElement.from_file(Z_BANDPASS_FILE, wave_unit='nm')

    #Define the source spectrum
    wave, flux, err = np.loadtxt(SPECTRA_FILE, unpack=True)
    wave = wave * u.angstrom
    flux = flux * (1e-17 * u.erg)/u.s/(u.cm**2)/u.angstrom
    spectrum = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)


    #Convolving the filters with the spectra
    narrow_observation = Observation(spectrum, narrow_band_pass)
    i_observation = Observation(spectrum, i_band_pass, force='taper')
    z_observation = Observation(spectrum, z_band_pass, force = 'taper')

    #Working out the AB magnitudes
    narrow_mag = narrow_observation.effstim(u.ABmag)
    i_mag = i_observation.effstim(u.ABmag)
    z_mag = z_observation.effstim(u.ABmag)

    #Working out the colors.
    i_n = i_mag - narrow_mag
    i_z = i_mag - z_mag
    n_z = narrow_mag - z_mag
