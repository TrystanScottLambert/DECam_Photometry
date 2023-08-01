"""
Script to plot the cummulative plot using the selection criteria from
Eduardo and Chira in 2013 and 2017.
"""
import numpy as np
import pylab as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from pandas import DataFrame

from plot_cumulative import cross_match_to_sexcat, add_ab_mags
from sex_catalog import SExtractorCat
from zero_points import zero_points, ZeroPoints
from zero_points_cdfs import zero_points_cdfs
from plotting import start_plot, end_plot


INNER_AREA_OF_SUPRESSION = 5487762.05698824
CDFS_AREA = 20385305 # Arcseconds
DECAM_AREA = 2.87e7 # Arcseconds 23212237.94301176  #The area minus the inner region.
IMACS_AREA = 5.01046e6 # Arcseconds

INFILE_US = 'candidates_e.txt'
SEX_US = '../correct_stacks/N964/n964.cat'
INFILE_CDFS = 'candidates_cdfs_e.txt'
SEX_CDFS = '../CDFS_LAGER/n964_cdfs.cat'

ra_us, dec_us = np.loadtxt(INFILE_US, unpack=True)
ra_cdfs, dec_cdfs = np.loadtxt(INFILE_CDFS, unpack=True)

our_cat = cross_match_to_sexcat(ra_us, dec_us, SExtractorCat(SEX_US))
cdfs_cat = cross_match_to_sexcat(ra_cdfs, dec_cdfs, SExtractorCat(SEX_CDFS))

our_cat['MAG_CORR'] = our_cat['MAG_APER'] + zero_points.n964_band.mag_correct(1)
cdfs_cat['MAG_CORR'] = cdfs_cat['MAG_APER'] + zero_points_cdfs.n964_band.mag_correct(1)

limits = np.arange(23., 25.1, 0.1)
our_counts = []
cdfs_counts = []
for limit in limits:
    our_counts.append(len(np.where(our_cat['MAG_CORR'] < limit)[0]))
    cdfs_counts.append(len(np.where(cdfs_cat['MAG_CORR'] < limit)[0]))

start_plot('magnitude limit', r'N(< mag) [arcseconds$^{-2}$]')
plt.plot(limits, np.array(our_counts) / DECAM_AREA, label = 'This Work')
plt.plot(limits, np.array(cdfs_counts) / CDFS_AREA, label = 'CDFS')
plt.legend()
end_plot('plots/cummulative_eduardo.png')

start_plot('magnitdue limit', 'Normalized-Count Ratio')
plt.plot(limits, (np.array(our_counts) / DECAM_AREA)/(np.array(cdfs_counts) / CDFS_AREA))
end_plot('plots/eduardo_ratio.png')

