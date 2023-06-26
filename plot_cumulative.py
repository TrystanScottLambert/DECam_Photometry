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


def cross_match_to_sexcat(
        ra_array: np.ndarray, dec_array:np.ndarray, sex_catalog: SExtractorCat
        ) -> DataFrame:
    """Reduces the sexcatalog to the cross matched values"""
    sex_ra, sex_dec = np.array(
        sex_catalog.catalog['ALPHAPEAK_J2000']), np.array(sex_catalog.catalog['DELTAPEAK_J2000'])
    catalog = SkyCoord(ra = sex_ra*u.deg, dec = sex_dec*u.deg)
    candidate_coords = SkyCoord(ra = ra_array * u.deg, dec = dec_array*u.deg)
    idx, _, _ = candidate_coords.match_to_catalog_sky(catalog)
    return sex_catalog.catalog.iloc[idx]


#n_cat = SExtractorCat('../correct_stacks/N964/n964.cat')
#n135_cat = SExtractorCat('../correct_stacks/N964/n964_135.cat')
n_cat = SExtractorCat('imacs_n964.cat')
n135_cat = SExtractorCat('imacs_n964_135.cat')

n_cdfs_cat = SExtractorCat('../CDFS_LAGER/n964_cdfs.cat')
n135_cdfs_cat = SExtractorCat('../CDFS_LAGER/n964_135_cdfs.cat')

candidates = 'candidates_imacs.txt'
#candidates = 'candidates.txt'
candidates_cdfs = 'candidates_cdfs.txt'

ra_us, dec_us = np.loadtxt(candidates, unpack=True)
ra_cdfs, dec_cdfs = np.loadtxt(candidates_cdfs, unpack=True)

matched_n_cat = cross_match_to_sexcat(ra_us, dec_us, n_cat)
matched_n135_cat = cross_match_to_sexcat(ra_us, dec_us, n135_cat)
matched_n_cat_cdfs = cross_match_to_sexcat(ra_cdfs, dec_cdfs, n_cdfs_cat)
matched_n135_cat_cdfs = cross_match_to_sexcat(ra_cdfs, dec_cdfs, n135_cdfs_cat)

### Narrowband
#ratio = 2.87e7/20385305
ratio = 5.01046e6/20385305
matched_n_cat['MAG_CORR'] = matched_n_cat['MAG_APER'] + zero_points.n964_band.mag_correct(1)
matched_n135_cat['MAG_CORR'] = matched_n135_cat['MAG_APER'] + zero_points.n964_band.mag_correct(1.35/2)
n135_lim = 25.10
n_lim = 24.66
depth = -np.arange(1, -0.01, -0.01)
counts = []
for i in depth:
    cut_n = np.where(matched_n_cat['MAG_CORR'] < n_lim + i)[0]
    cut_135 = np.where(matched_n135_cat['MAG_CORR'] < n135_lim + i)[0]
    counts.append(len(np.intersect1d(cut_n, cut_135)))

matched_n_cat_cdfs['MAG_CORR'] = matched_n_cat_cdfs['MAG_APER'] + zero_points_cdfs.n964_band.mag_correct(1)
matched_n135_cat_cdfs['MAG_CORR'] = matched_n135_cat_cdfs['MAG_APER'] + zero_points_cdfs.n964_band.mag_correct(1.35/2)
counts_cdfs = []
for i in depth:
    cut_n = np.where(matched_n_cat_cdfs['MAG_CORR'] < n_lim + i)[0]
    cut_135 = np.where(matched_n135_cat_cdfs['MAG_CORR'] < n135_lim + i)[0]
    counts_cdfs.append(len(np.intersect1d(cut_n, cut_135)))

start_plot('Relative Narrow Band Selection', r'Counts per Area [arcsecond$^{-2}$]')
plt.plot(depth, np.array(counts)/5.01046e6, label='This work', color='k', lw=3)
plt.plot(depth, np.array(counts_cdfs)/20385305, label='CDFS', color='r', ls=':', lw=2)
plt.legend()
end_plot('cumulative_plot.png')
