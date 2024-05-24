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
from synphot import SourceSpectrum, etau_madau, ReddeningLaw
from synphot.models import Empirical1D
import synphot.units as su


cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def gaussian(x_array: np.ndarray, mean: float, sigma: float, amplitude: float) -> np.ndarray:
    """Definition of the gaussian function"""
    return amplitude * np.exp(-np.power(x_array - mean, 2.) / (2 * np.power(sigma, 2.)))

#try 1
class SedSpectrum:
    """SED spectrum provided by https://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html"""

    def __init__(self, file_name: str, redshift: float, extinction: float = 0) -> None:
        self.rest_spectrum = SourceSpectrum.from_file(
            file_name, keep_neg=True, wave_unit='Angstrom', flux_unit='FLAM')
        self.red_shifted_spectrum = SourceSpectrum(
            self.rest_spectrum.model, z=redshift, z_type='conserve_flux')
        self.wave = self.rest_spectrum.to_spectrum1d().spectral_axis
        self.extinction = extinction
        self.z = redshift
        self.flux = self.rest_spectrum.to_spectrum1d().flux

    @property
    def spectrum(self):
        """Redshifted, and if appropriate, extinction corrected spectrum."""
        flux = self.flux/(1+self.z)
        spectrum = SourceSpectrum(Empirical1D, points=self.wave*(1+self.z), lookup_table=flux, keep_neg=True)
        extinction_curve = ReddeningLaw.from_extinction_model('xgalsb').extinction_curve(self.extinction)
        spectrum = spectrum * extinction_curve
        return spectrum

class LaeSpectrum:
    """Class representation of a lyman alpha emitter spectrum."""
    LLYA_REST = 1215.67 * u.angstrom
    FWHM_LYA  = 200.    * (u.km/u.s)
    SIG_LYA   = (LLYA_REST * (FWHM_LYA/c.c.to(u.km/u.s))) * gaussian_fwhm_to_sigma # Ang.
    EW_LYA    = 50. *u.angstrom
    LOG_LUM_LYA = 43.33 #Matthee+15, Table 6, alpha=-1.5 UDS+COSMOS+SA22 (flux in erg/s)
    wavelength_rest = np.arange(500, 2500, 1) * u.angstrom
    flux_rest = gaussian(wavelength_rest, LLYA_REST, SIG_LYA, 1.)

    def __init__(self, redshift: float, include_absorption=True) -> None:
        """Initializaing the spectrum based on the redshift"""
        self.lum_dist = cosmo.luminosity_distance(redshift).to(u.cm)
        self._log_flux = LaeSpectrum.LOG_LUM_LYA - np.log10(
            4 * np.pi * self.lum_dist.value**2) #- 48
        self._flux = 10**self._log_flux * (u.erg/u.s/u.cm**2)
        self._ew_lya_obs = LaeSpectrum.EW_LYA * (1 + redshift)

        wavelength = LaeSpectrum.wavelength_rest * (1 + redshift)
        flux = (LaeSpectrum.flux_rest.value * self._flux.value) + (self._flux.value/self._ew_lya_obs.value)
        #flux = (LaeSpectrum.flux_rest * self._flux/1) + (self._flux/self._ew_lya_obs)
        #flux = su.convert_flux(wavelength, flux, out_flux_unit=su.PHOTLAM)

        self.spectrum = SourceSpectrum(
            Empirical1D, points = wavelength, lookup_table = flux * (su.FLAM), keep_neg=True)
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

if __name__ == '__main__':
    qso_spectrum = LaeSpectrum(6.9)
    qso_spectrum.spectrum.plot()
    plt.show()
    