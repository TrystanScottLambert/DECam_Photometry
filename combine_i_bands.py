"""
Combining the imacs data and decam data for the i band 
using an average weighting, in order to determine
a deeper magnitude limit.
"""

import numpy as np
import pylab as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from photutils.aperture import SkyCircularAperture, ApertureStats


from sex_catalog import SExtractorCat
from zero_points import zero_points

IMACS_ZPT = 27.633

def get_imacs_datawcs(imacs_name: str) -> fits.HDUList:
    """Reads in the imacs data. Returns the data and wcs"""
    hdul = fits.open(imacs_name)
    return hdul[0].data, WCS(hdul[0].header)

def perform_imacs_photometry(
        imacs_data: np.ndarray, imacs_wcs: WCS,
        apertures: SkyCircularAperture, uncertainty_map: np.ndarray
        ) -> np.ndarray:
    """Does aperture photometry, counts to instrumental mags, to AB mags, to flux density"""
    phot_table = ApertureStats(imacs_data, apertures, wcs = imacs_wcs, error=uncertainty_map)
    counts = phot_table.sum
    count_err = phot_table.sum_err
    mag_err = propagate_conterr_magerr(counts, count_err)
    flux_density = convert_counts_to_flux_density(counts, IMACS_ZPT)
    flux_density_err = propgate_magerr_fluxerr(flux_density, mag_err)
    return flux_density, flux_density_err


def convert_abmag_to_flux_density(ab_mags: np.ndarray) -> np.ndarray:
    """Converts the abmags to flux density."""
    return 10**((ab_mags + 48.6)/(-2.5)) * u.erg * (u.cm**(-2)) * (u.s**(-1)) * (u.Hz**(-1))

def propagate_conterr_magerr(counts: np.ndarray, count_errors: np.ndarray) -> np.ndarray:
    """
    Converts the count erros into instrumental magnitude errors.
    see: https://sextractor.readthedocs.io/en/latest/Param.html#fluxes-and-magnitudes
    """
    return (2.5/np.log(10)) * (count_errors/counts)

def propgate_magerr_fluxerr(flux_density: np.ndarray, mag_uncertainty: np.ndarray) -> np.ndarray:
    """Converts the uncertatinty of the magnitudes into the uncertainty for flux."""
    return flux_density * np.abs((np.log(10) * mag_uncertainty)/2.5)

def convert_counts_to_flux_density(counts: np.ndarray, zero_point: float):
    """Converts counts into flux density."""
    instrumental_mags = -2.5*np.log10(counts)
    mags = instrumental_mags + zero_point
    flux_density = convert_abmag_to_flux_density(mags)
    return flux_density

def calculate_weighted_flux(fluxes: np.ndarray, flux_uncertainties: np.ndarray) -> np.ndarray:
    """
    Determines the weighted flux for the given bands.
    Arrays represent the flux of a single source in 
    multiple bands and not multiple sources in a single band.
    """
    return np.sum(fluxes/(flux_uncertainties**2)) / np.sum(1./(flux_uncertainties**2))

def calculate_weighted_fluxerr(flux_uncertainties: np.ndarray) -> np.ndarray:
    """Determines the flux uncertainties of a single source across multiple bands."""
    return np.sum((1./(flux_uncertainties**2))**(-0.5))

if __name__ == '__main__':
    DECAM_CAT_NAME  = '../correct_stacks/N964/i.cat'
    IMACS_NIGHT1 = '../../IMACS_photometry/imacs_data/night_1_theli.fits'
    IMACS_NIGHT1_WEIGHT = '../../IMACS_photometry/imacs_data/night_1_theli.weight.fits'
    IMACS_NIGHT2 = '../../IMACS_photometry/imacs_data/night_2_theli.fits'
    IMACS_NIGHT2_WEIGHT = '../../IMACS_photometry/imacs_data/night_2_theli.weight.fits'

    hdu1_data, hdu1_wcs = get_imacs_datawcs(IMACS_NIGHT1)
    hdu2_data, hdu2_wcs = get_imacs_datawcs(IMACS_NIGHT2)
    hdu1_weight, _ = get_imacs_datawcs(IMACS_NIGHT1_WEIGHT)
    hdu2_weight, _ = get_imacs_datawcs(IMACS_NIGHT2_WEIGHT)

    uncertainty_1 = np.sqrt(1./hdu1_weight)
    uncertainty_2 = np.sqrt(1./hdu2_weight)

    #prepare decam data
    decam_catalog = SExtractorCat(DECAM_CAT_NAME)
    decam_ra = list(decam_catalog.catalog['ALPHAPEAK_J2000'])
    decam_dec = list(decam_catalog.catalog['DELTAPEAK_J2000'])
    decam_mags = np.array(decam_catalog.catalog['MAG_APER']) + zero_points.i_band.mag_correct(1)
    decam_mag_err = np.array(decam_catalog.catalog['MAGERR_APER'])

    decam_flux_density = convert_abmag_to_flux_density(decam_mags)
    decam_flux_density_err = propgate_magerr_fluxerr(decam_flux_density, decam_mag_err)
    positions = SkyCoord(ra = decam_ra * u.deg, dec = decam_dec * u.deg)
    apertures = SkyCircularAperture(positions, r=1. * u.arcsec)

    night1_flux_density = perform_imacs_photometry(hdu1_data, hdu1_wcs, apertures, uncertainty_1)
    night2_flux_density = perform_imacs_photometry(hdu2_data, hdu2_wcs, apertures, uncertainty_2)

    fluxes = np.dstack([night1_flux_density[0], night2_flux_density[0], decam_flux_density])[0]
    flux_errs = np.dstack([night1_flux_density[1], night2_flux_density[1], decam_flux_density_err])[0]

    weighted_fluxes = [calculate_weighted_flux(flux, flux_errs[i]) for i, flux in enumerate(fluxes)]
    weighted_flux_errs = [calculate_weighted_fluxerr(fluxerrs) for fluxerrs in flux_errs]
    F_NU = 3.631e-20 * u.erg * (u.cm**(-2)) * (u.s**(-1)) * (u.Hz**(-1)) # See https://en.wikipedia.org/wiki/AB_magnitude
    weighted_mags = np.array([-2.5 * np.log10(flux/F_NU) for flux in weighted_fluxes])

    snr = np.array([weighted_flux/weighted_flux_errs[i] for i, weighted_flux in enumerate(weighted_fluxes)])

    cut = np.where((snr>2.9) & (snr<3.1))[0]
    plt.hist(weighted_mags[cut])
    plt.axvline(np.median(weighted_mags[cut]),ls='--', color='k')
    plt.axvline(np.mean(weighted_mags[cut]))
    print('Depth value in DECAm mags is: ', np.median(weighted_mags[cut]))
    plt.xlabel('Weighted Magnitudes')
    plt.ylabel('Counts')
    plt.show()

    plt.scatter(weighted_mags, snr, s=0.1, color='k')
    plt.axhline(3, ls=':')
    plt.ylim(0, 15)
    plt.xlabel('Weighted Magnitudes')
    plt.ylabel('SNR')
    plt.show()