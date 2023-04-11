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

def get_matches(target_catalog: SkyCoord, main_catalog: SkyCoord):
    """
    Finds the sources in the target catalog which have appropriate matches
    in the main_catalog. This is determined through a tolerance of 1".

    Returns the indicies of the main_catalog and then the indicies of the 
    target catalog.
    """
    idx, d2d, _ = target_catalog.match_to_catalog_sky(main_catalog)
    cut = np.where(d2d < 3600 *u.arcsec)[0]
    plt.hist(d2d[cut])
    plt.show()
    return idx[cut], cut


if __name__ == '__main__':
    MOCK_SOURCES_2_SEXCAT = '../correct_stacks/N964/n964_false.cat'
    MOCK_SOURCES_135_SEXCAT = '../correct_stacks/N964/n964_false_135.cat'
    ORIGINAL_SEXCAT = '../correct_stacks/N964/n964.cat'
    MOCK_SOURCES_FILE = 'mock_lae_sources.txt'
    FITS_FILE = '../correct_stacks/N964/n964.injected.fits'

    with open(ORIGINAL_SEXCAT, encoding='utf8') as file:
        original_lines = file.readlines()

    with open(MOCK_SOURCES_2_SEXCAT, encoding='utf8') as file:
        mock_lines = file.readlines()


    ra_original, dec_original, flux = np.loadtxt(ORIGINAL_SEXCAT, usecols=(0, 1, 6), unpack=True)
    ra_mock, dec_mock, flux_mock = np.loadtxt(MOCK_SOURCES_2_SEXCAT, usecols=(0,1,6), unpack=True)

    original_cat = SkyCoord(ra = ra_original*u.deg, dec = dec_original*u.deg)
    mock_cat = SkyCoord(ra = ra_mock*u.deg, dec = dec_mock*u.deg)
    idx, d2d, _ = original_cat.match_to_catalog_sky(mock_cat)
