"""
Module for modelling and creating LAE spectra.
Creates a synphot spectrum of a LAE for a given redshift.
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy import constants as c
from astropy.cosmology import FlatLambdaCDM
from astropy.stats import gaussian_fwhm_to_sigma
from synphot import SourceSpectrum, etau_madau
from synphot.models import Empirical1D
import synphot.units as su


cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def gaussian(x_array: np.ndarray, mean: float, sigma: float, amplitude: float) -> np.ndarray:
    """Definition of the gaussian function"""
    return amplitude * np.exp(-np.power(x_array - mean, 2.) / (2 * np.power(sigma, 2.)))


class LaeSpectrum:
    """Class representation of a lyman alphe emitter spectrum."""
    LLYA_REST = 1215.67 * u.angstrom
    FWHM_LYA  = 200.    * (u.km/u.s)
    SIG_LYA   = (LLYA_REST * (FWHM_LYA/c.c.to(u.km/u.s))) * gaussian_fwhm_to_sigma # Ang.
    EW_LYA    = 50. *u.angstrom
    LOG_LUM_LYA = 43.33 #Matthee+15, Table 6, alpha=-1.5 UDS+COSMOS+SA22 (flux in erg/s)
    wavelength_rest = np.arange(500, 2500, 1) * u.angstrom
    flux_rest = gaussian(wavelength_rest, LLYA_REST, SIG_LYA, 1./u.angstrom)

    def __init__(self, redshift: float, include_absorption=True) -> None:
        """Initializaing the spectrum based on the redshift"""
        self.lum_dist = cosmo.luminosity_distance(redshift).to(u.cm)
        self._log_flux = LaeSpectrum.LOG_LUM_LYA - np.log10(
            4 * np.pi * self.lum_dist.value**2) - 48
        self._flux = 10**self._log_flux * (u.erg/u.s/self.lum_dist.unit**2)


        self._ew_lya_obs = LaeSpectrum.EW_LYA * (1 + redshift)

        wavelength = LaeSpectrum.wavelength_rest * (1 + redshift)
        flux = (LaeSpectrum.flux_rest * self._flux/1) + (self._flux/self._ew_lya_obs)
        flux = su.convert_flux(wavelength, flux, out_flux_unit=su.PHOTLAM)


        self.spectrum = SourceSpectrum(
            Empirical1D, points = wavelength, lookup_table = flux, keep_neg=True)
        if include_absorption:
            extinction_curve = etau_madau(wavelength, redshift)
            self.spectrum = self.spectrum * extinction_curve

    def save_plot(self, outfile: str) -> None:
        """Saves a figure of the plot."""
        self.spectrum.plot()
        plt.xlim(7000, 11500)
        plt.ylim([0.,(self._flux/self._ew_lya_obs)*10.9])
        plt.savefig(outfile)

    def to_file(self, outfile: str) -> None:
        """Creates a text file of the spectrum."""
        spectrum_1d = self.spectrum.to_spectrum1d()
        wavelengths = spectrum_1d.spectral_axis.value
        fluxes = spectrum_1d.flux.value
        with open(outfile, 'w', encoding='utf8') as output_file:
            output_file.write(f'#{spectrum_1d.spectral_axis.unit} {spectrum_1d.flux.unit} \n')
            for wavelength, flux in zip(wavelengths, fluxes):
                try:
                    output_file.write(f'{wavelength} {flux} \n')
                except IndexError:
                    pass
