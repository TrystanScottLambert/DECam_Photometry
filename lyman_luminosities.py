"""This module is for determining the lyman-alpha luminosity from NB964 filter and z-filter."""

import numpy as np
import pylab as plt
from synphot import etau_madau
from synphot import SpectralElement
import astropy.constants as cons
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
import astropy.units as u

from zero_points_cdfs import zero_points_cdfs

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
        return self.transmission.avgwave()

    @property
    def area_of_filter(self):
        """Determines the integral of the curve. Since these are transmission curves the integration 
        will be in units of angstrom (as apposed to angstrom Flux for example.)"""
        return self.transmission.integrate()

    @property
    def filter_constant(self):
        """Determines the filter constant from Eq 4, Lambert et.al (2023)"""
        return cons.c * (self.area_of_filter / (self.central_wavelength**2))

    @property
    def transmission_at_lyman(self):
        """The transmission at the observed wavelength of lyman alpha."""
        return self.transmission(LYMAN_ALPHA_OBSERVED_WAVELENGTH)

    @property
    def continuum_integral(self):
        """The integral seen in equation 7 of lambert et. al."""
        wave = self.transmission.to_spectrum1d().wavelength
        attenuated_transmission = etau_madau(wave, QSO_REDSHIFT) * self.transmission
        attenuated_values = attenuated_transmission.to_spectrum1d().flux
        integrand = attenuated_values * (wave**(-2))
        return numerical_integration(wave.value, integrand.value)*(wave.unit**-1)

def calculate_lyman_alpha_flux(
        measured_nb964_flux: float, measured_z_flux: float, nb964: Filter, z: Filter,
        measured_nb964_err: float, measured_z_flux_err: float):
    """
    Determines the lyman alpha flux from the narrowband and broadband measured flux values.
    see http://physics.uwyo.edu/~chip/Classes/ASTR4610/Lec_Distances.pdf for units.
    """
    ratio = z.continuum_integral / nb964.continuum_integral
    numerator_t1 = z.filter_constant * measured_z_flux
    numerator_t2 = nb964.filter_constant * measured_nb964_flux * ratio
    denominator_t1 = z.transmission_at_lyman
    denominator_t2 = nb964.transmission_at_lyman * ratio
    value = (numerator_t1 - numerator_t2)/(denominator_t1 - denominator_t2)
    
    #uncertainties
    u_n1 = z.filter_constant * measured_z_flux_err
    u_n2 = nb964.filter_constant*measured_nb964_err*ratio
    uncertainty = np.hypot(u_n1, u_n2)/(denominator_t1-denominator_t2)
    return value.to(u.erg * (u.s**(-1)) * (u.cm**(-2))), np.abs(uncertainty.to(u.erg * (u.s**(-1)) * (u.cm**(-2))))

def convert_flux_to_luminosity(flux: float):
    """converts flux into luminosity"""
    lum_distance = COSMO.luminosity_distance(QSO_REDSHIFT)
    value = (flux * np.pi * 4 * (lum_distance**2))#/(QSO_REDSHIFT + 1)
    return value.to(u.erg/u.s) #converting to same units as hu et. al., 2019


#preparing the catalog information.
def write_file(
        outfile_name: str, lya_lums: np.ndarray, lya_err: np.ndarray,
        sfrs: np.ndarray, sfrs_err: np.ndarray,
        mag: np.ndarray, mag_err: np.ndarray) -> None:
    """Writing the calculated luminosities and star formation rates to file."""
    with open(outfile_name, 'w', encoding='utf8') as file:
        for lya, l_err, sfr, s_err, m, me in zip(lya_lums, lya_err, sfrs, sfrs_err, mag, mag_err):
            file.write(f'{m} {me} {lya.value/1e42} {l_err.value/1e42} {sfr} {s_err} \n')


if __name__ == '__main__':
    NARROW_BANDPASS_FILE = '../QSO_Spectra/NB964_DECam_F29.txt'
    Z_BANDPASS_FILE = '../QSO_Spectra/decam_z_bandpass.txt'

    NARROW_CATALOG = '../correct_stacks/N964/n964.cat'
    Z_CATALOG = '../correct_stacks/N964/z.cat'
    CANDIDATES_CATALOG = 'candidates_e.txt'
    OUTFILE = 'candidate_luminosities.dat'
    
    #NARROW_CATALOG = '../CDFS_LAGER/n964_cdfs.cat'
    #Z_CATALOG = '../CDFS_LAGER/z_cdfs.cat'
    #CANDIDATES_CATALOG = 'candidates_cdfs_e.txt'
    #OUTFILE = 'candidate_luminosities_cdfs.dat'


    ra_n, dec_n, n_mag, n_err = np.loadtxt(NARROW_CATALOG, unpack=True, usecols=(0, 1, 4, 5))
    z_mag, z_err = np.loadtxt(NARROW_CATALOG, unpack=True, usecols=(4,5))

    ra_candidates, dec_candidates = np.loadtxt(CANDIDATES_CATALOG, unpack=True)
    candidates = SkyCoord(ra = ra_candidates * u.deg, dec = dec_candidates * u.deg)
    n_catalog = SkyCoord(ra = ra_n * u.deg, dec = dec_n * u.deg)
    idx, d2d, _ = candidates.match_to_catalog_sky(n_catalog)

    n_mag = n_mag[idx] + zero_points_cdfs.n964_band.mag_correct(1)
    z_mag = z_mag[idx] + zero_points_cdfs.z_band.mag_correct(1)
    n_err = n_err[idx]
    z_err = z_err[idx]

    '''
    See table at https://pysynphot.readthedocs.io/en/latest/units.html 
    for units of fnu which is what is measured photometrically.
    '''
    n_flux_nu = 10**((n_mag + 48.6)/(-2.5)) * u.erg * (u.cm**(-2)) * (u.s**(-1)) * (u.Hz**(-1))
    z_flux_nu = 10**((z_mag + 48.6)/(-2.5)) * u.erg * (u.cm**(-2)) * (u.s**(-1)) * (u.Hz**(-1))
    n_flux_nu_err = np.log(10)*n_flux_nu*(n_err/2.5)
    z_flux_nu_err = np.log(10)*z_flux_nu*(z_err/2.5)
    

    NB964 = Filter(NARROW_BANDPASS_FILE)
    Z_BAND = Filter(Z_BANDPASS_FILE)

    lya_flux, lya_flux_uncertainty = calculate_lyman_alpha_flux(n_flux_nu, z_flux_nu, NB964, Z_BAND, n_flux_nu_err, z_flux_nu_err)
    lya_lum = convert_flux_to_luminosity(lya_flux)
    lya_lum_err = convert_flux_to_luminosity(lya_flux_uncertainty) # assuming no error on luminosity distance
    log_10_lya = np.log10(lya_lum.value)

    ha_lum = lya_lum / 8.7 # See muzzucchelli 2017
    ha_lum_err = lya_lum_err/8.7
    log_10_sfr = np.log10(ha_lum.value) - 41.27 #See muzzucchelli 2017
    log_10_sfr_err = (1./np.log(10)) * (ha_lum_err/ha_lum)

    sfr = 10**(log_10_sfr)
    sfr_err = np.log(10)*sfr*log_10_sfr_err

    write_file(OUTFILE, lya_lum, lya_lum_err, sfr, sfr_err, n_mag, n_err)
    plt.hist(log_10_lya, bins=np.arange(42.0, 43.8,0.1))
    plt.show()
