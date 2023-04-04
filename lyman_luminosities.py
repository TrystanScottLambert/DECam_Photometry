"""This module is for determining the lyman-alpha luminosity from NB964 filter and z-filter."""

import numpy as np
from synphot import etau_madau
from synphot import SpectralElement
import astropy.constants as cons
from astropy.cosmology import FlatLambdaCDM

COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
QSO_REDSHIFT = 6.9018
LYMAN_ALPHA_OBSERVED_WAVELENGTH = 1216 * (QSO_REDSHIFT + 1) # angstrom


def numerical_integration(x_array, y_array):
    """
    Simple approximation to numerical integrate over arrays (instead of over functions).
    These will be more accurate the smaller the differences in x[i] and x[i+1].
    """
    x_bins = np.array([x_array[i+1] - x_array[i] for i in range(len(x_array)-1)])
    y_avg = np.array([(y_array[i+1] + y_array[i])/2 for i in range(len(y_array)-1)])
    return np.sum(x_bins * y_avg)


class Filter:
    """Properties of the filters"""
    def __init__(self, tranmission_file: str):
        self.transmission = SpectralElement.from_file(tranmission_file, wave_unit='nm')

    @property
    def central_wavelength(self):
        """Determines the average wavelength of the filter in angstroms"""
        return self.transmission.avgwave().value

    @property
    def area_of_filter(self):
        """Determines the integral of the curve. Since these are transmission curves the integration 
        will be in units of angstrom (as apposed to angstrom Flux for example.)"""
        return self.transmission.integrate().value

    @property
    def filter_constant(self):
        """Determines the filter constant from Eq 4, Lambert et.al (2023)"""
        return cons.c.value * (self.area_of_filter / (self.central_wavelength**2))

    @property
    def transmission_at_lyman(self):
        """The transmission value at the observed wavelength of lyman alpha."""
        return self.transmission(9608).value

    @property
    def continuum_integral(self):
        """The integral seen in equation 7 of lambert et. al."""
        wave = self.transmission.to_spectrum1d().wavelength.value
        attenuated_transmission = etau_madau(wave, QSO_REDSHIFT) * self.transmission
        attenuated_values = attenuated_transmission.to_spectrum1d().flux.value
        integrand = attenuated_values * (wave**(-2))
        return numerical_integration(wave, integrand)

def calculate_c(measured_nb964_flux: float, measured_z_flux: float, nb964: Filter, z: Filter):
    """Determines the C constant in the similtaneous equations"""
    term_1 = z.filter_constant * measured_z_flux * nb964.transmission_at_lyman 
    term_2 = nb964.filter_constant * measured_nb964_flux * z.transmission_at_lyman
    return (term_1 - term_2) / z.continuum_integral

def calculate_lyman_alpha_flux(measured_nb964_flux: float, measured_z_flux: float, nb964: Filter, z: Filter):
    """Determines the lyman alpha flux from the narrowband and broadband measured flux values."""
    constant = calculate_c(measured_nb964_flux, measured_z_flux, nb964, z)
    term_1 = nb964.filter_constant * measured_nb964_flux
    term_2 = constant * nb964.continuum_integral
    return (term_1 - term_2)/nb964.transmission_at_lyman

def convert_flux_to_luminosity(flux: float):
    """converts flux into luminosity"""
    lum_distance = COSMO.luminosity_distance(QSO_REDSHIFT)
    return flux * np.pi * 4 * (lum_distance**2)

if __name__ == '__main__':
    SPECTRA_FILE = '../QSO_Spectra/J2348_dered.txt'
    NARROW_BANDPASS_FILE = '../QSO_Spectra/NB964_DECam_F29.txt'
    I_BANDPASS_FILE = '../QSO_Spectra/decam_i_bandpass.txt'
    Z_BANDPASS_FILE = '../QSO_Spectra/decam_z_bandpass.txt'


    NB964 = Filter(NARROW_BANDPASS_FILE)

    #Define the filter objects
    narrow_band_pass = SpectralElement.from_file(NARROW_BANDPASS_FILE, wave_unit='nm')
    i_band_pass = SpectralElement.from_file(I_BANDPASS_FILE, wave_unit='nm')
    z_band_pass = SpectralElement.from_file(Z_BANDPASS_FILE, wave_unit='nm')
