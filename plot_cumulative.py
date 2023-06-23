"""
Making the cummulative plot to show the overdensity.
"""

import numpy as np
import pylab as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from pandas import DataFrame

from sex_catalog import SExtractorCat
from zero_points import zero_points
from zero_points_cdfs import zero_points_cdfs
from plotting import start_plot, end_plot


def cross_match_to_sexcat(ra_array: np.ndarray, dec_array:np.ndarray, sex_catalog: SExtractorCat) -> DataFrame:
    """Reduces the sexcatalog to the cross matched values"""
    sex_ra, sex_dec = np.array(sex_catalog.catalog['ALPHAPEAK_J2000']), np.array(sex_catalog.catalog['DELTAPEAK_J2000'])
    catalog = SkyCoord(ra = sex_ra*u.deg, dec = sex_dec*u.deg)
    candidate_coords = SkyCoord(ra = ra_array * u.deg, dec = dec_array*u.deg)
    idx, _, _ = candidate_coords.match_to_catalog_sky(catalog)
    return sex_catalog.catalog.iloc[idx]


i_cat = SExtractorCat('../correct_stacks/N964/i.cat')
z_cat = SExtractorCat('../correct_stacks/N964/z.cat')
n_cat = SExtractorCat('../correct_stacks/N964/n964.cat')
n135_cat = SExtractorCat('../correct_stacks/N964/n964_135.cat')

i_cdfs_cat = SExtractorCat('../CDFS_LAGER/i_cdfs.cat')
z_cdfs_cat = SExtractorCat('../CDFS_LAGER/z_cdfs.cat')
n_cdfs_cat = SExtractorCat('../CDFS_LAGER/n964_cdfs.cat')
n135_cdfs_cat = SExtractorCat('../CDFS_LAGER/n964_135_cdfs.cat')

candidates = 'candidates.txt'
candidates_cdfs = 'candidates_cdfs.txt'

ra_us, dec_us = np.loadtxt(candidates, unpack=True)
ra_cdfs, dec_cdfs = np.loadtxt(candidates_cdfs, unpack=True)

matched_i_cat = cross_match_to_sexcat(ra_us, dec_us, i_cat)
matched_z_cat = cross_match_to_sexcat(ra_us, dec_us, z_cat)
matched_n_cat = cross_match_to_sexcat(ra_us, dec_us, n_cat)
matched_n135_cat = cross_match_to_sexcat(ra_us, dec_us, n135_cat)

matched_i_cat_cdfs = cross_match_to_sexcat(ra_cdfs, dec_cdfs, i_cdfs_cat)
matched_z_cat_cdfs = cross_match_to_sexcat(ra_cdfs, dec_cdfs, z_cdfs_cat)
matched_n_cat_cdfs = cross_match_to_sexcat(ra_cdfs, dec_cdfs, n_cdfs_cat)
matched_n135_cat_cdfs = cross_match_to_sexcat(ra_cdfs, dec_cdfs, n135_cdfs_cat)

### Manual method:
ratio = 2.87e7/20385305
x = -np.arange(1, -0.1, -0.1)
us = np.array([18, 22, 25, 29, 38, 56, 82, 121, 186, 251, 310])
them = np.array([16, 17, 20, 22, 28, 35, 55, 82, 137, 222, 309])
plt.plot(x, us, label='this work', lw=3, color='r')
plt.plot(x, them * ratio, label='Hu', lw=2, ls=':', color='k')
plt.legend()
plt.show()