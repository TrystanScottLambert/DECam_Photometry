"""
Module which determines the completeness value based on mock source injection. 
See Figure 4, section 4.1, and f_{comp} of equation 5 in Hu et. al., (2019).
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import pylab as plt
from inject_false_sources import MIN_MAG, MAX_MAG
from hu_2019_plot import FILTERS

def get_matches(target_catalog: SkyCoord, main_catalog: SkyCoord):
    """
    Finds the sources in the target catalog which have appropriate matches
    in the main_catalog. This is determined through a tolerance of 1".

    Returns the indicies of the main_catalog and then the indicies of the 
    target catalog.
    """
    idx, d2d, _ = target_catalog.match_to_catalog_sky(main_catalog)
    cut = np.where(d2d < 1 *u.arcsec)[0]
    return cut, idx[cut]


def remove_original_sources(recovered_catalog: SkyCoord, original_catalog: SkyCoord):
    """
    Matches the sources that are identical in the original catalog and removes them.
    If the positions are exactly the same as those in the original catalog then they have 
    to be sources from the original catalog.
    """

    _, d2d,_ = recovered_catalog.match_to_catalog_sky(original_catalog)
    cut = np.where(d2d != 0)[0]
    return cut

if __name__ == '__main__':
    MOCK_SOURCES_2_SEXCAT = '../correct_stacks/N964/n964_false.cat'
    MOCK_SOURCES_135_SEXCAT = '../correct_stacks/N964/n964_false_135.cat'
    ORIGINAL_SEXCAT = '../correct_stacks/N964/n964.cat'
    MOCK_SOURCES_FILE = 'mock_lae_sources.txt'
    FITS_FILE = '../correct_stacks/N964/n964.injected.fits'

    hdu = fits.open(FITS_FILE)
    wcs = WCS(hdu[0].header)
    x, y, mag = np.loadtxt(MOCK_SOURCES_FILE, unpack=True)
    flux = 10**((mag - FILTERS['n964'].zpt)/(-2.5))
    ra, dec = wcs.pixel_to_world_values(x, y)

    ra_original, dec_original, flux_original = np.loadtxt(ORIGINAL_SEXCAT, usecols=(0, 1, 6), unpack=True)
    ra_mock, dec_mock, x_mock, y_mock, mag_mock, flux_mock = np.loadtxt(MOCK_SOURCES_2_SEXCAT, usecols=(0, 1, 2, 3, 4, 6), unpack=True)
    mag_mock += FILTERS['n964'].zpt

    original_cat = SkyCoord(ra = ra_original * u.deg, dec = dec_original * u.deg)
    recovered_cat = SkyCoord(ra = ra_mock * u.deg, dec = dec_mock * u.deg)
    cat = SkyCoord(ra = ra * u.deg, dec = dec * u.deg)

    # Match the injected catalog with the known mock positions.
    idx_recovered, idx_mock = get_matches(recovered_cat, cat)

    # Remove any sources that are actually from the original catalog
    possible_recovered_idx = remove_original_sources(recovered_cat[idx_recovered], original_cat)


    plt.scatter(flux[idx_mock][possible_recovered_idx], flux_mock[idx_recovered][possible_recovered_idx])
    plt.plot(np.arange(min(flux), max(flux)), np.arange(min(flux), max(flux)), color='k', ls='--', lw=3)
    plt.xlabel('Mock Sources')
    plt.ylabel('SExtractor recovered Sources')
    plt.show()


    BINWIDTH = 0.2
    bins = np.arange(MIN_MAG, MAX_MAG + BINWIDTH, BINWIDTH)
    x_avg = [(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]

    y_mock, _ = np.histogram(mag, bins=bins)
    y_recovered, _ = np.histogram(mag[idx_mock], bins=bins)

    y_completeness = y_recovered/y_mock
    plt.plot(x_avg, y_completeness)
    plt.show()
