"""
Script to plot the cummulative plot using the selection criteria from
Eduardo and Chira in 2013 and 2017.
"""
import numpy as np
import pylab as plt
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from pandas import DataFrame
from astropy.io import fits

from plot_cumulative import cross_match_to_sexcat, add_ab_mags
from sex_catalog import SExtractorCat
from zero_points import zero_points, ZeroPoints
from zero_points_cdfs import zero_points_cdfs
from plotting import start_plot, end_plot


RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24) * u.deg
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1 *u.deg
qso_position = SkyCoord(ra = RA_QSO, dec = DEC_QSO)
cdfs_hdu = fits.open('CDFS_MASK.fits')
pix_scale = cdfs_hdu[0].header['PC2_2'] * 3600 #arcseconds


INNER_AREA_OF_SUPRESSION = 5487762.05698824
CDFS_AREA = len(np.where(cdfs_hdu[0].data == 1)[0]) * (pix_scale**2) # Arcseconds
DECAM_AREA = 2.87e7 # Arcseconds 23212237.94301176  #The area minus the inner region.
IMACS_AREA = 5.01046e6 # Arcseconds

INFILE_US = 'candidates_e.txt'
SEX_US = '../correct_stacks/N964/n964.cat'
INFILE_CDFS = 'candidates_cdfs_e.txt'
SEX_CDFS = '../CDFS_LAGER/n964_cdfs.cat'

ra_us, dec_us = np.loadtxt(INFILE_US, unpack=True)
c = SkyCoord(ra=ra_us*u.deg, dec = dec_us*u.deg)
distances_to_quasar = c.separation(qso_position)
ra_cdfs, dec_cdfs = np.loadtxt(INFILE_CDFS, unpack=True)

REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
ARCSEC_PER_KPC = COSMO.arcsec_per_kpc_comoving(REDSHIFT_QSO)
DEG_PER_MPC = ARCSEC_PER_KPC.to(u.deg / u.Mpc)
inner_region_distance = DEG_PER_MPC * 75 * u.Mpc
results = distances_to_quasar < inner_region_distance
inner = len(np.where(results ==True)[0])
outer = len(np.where(results == False)[0])
inner_area = (inner_region_distance.to(u.arcsec) **2) * np.pi
outer_area = 23212237.94301176




#Printing values
print('CDFS density [arcseconds ^{-2}]: ', len(ra_cdfs)/CDFS_AREA)
print('our density [arcseconds ^{-2})]: ', len(ra_us)/DECAM_AREA)
print('Ratio: ', (len(ra_us)/DECAM_AREA)/(len(ra_cdfs)/CDFS_AREA))
print('-----------------------------')
print('Inner 75cMpc is:')
print('Density[arcseconds ^{-2})]: ', inner/inner_area.value)
print('Ratio: ', (inner/inner_area.value)/(len(ra_cdfs)/CDFS_AREA))
print('-----------------------------')
print('Outer 75 cMpc')
print('Density[arcseconds ^{-2})]: ', outer/23212237.94301176)
print('Ratio: ', (outer/23212237.94301176)/(len(ra_cdfs)/CDFS_AREA))
print('-----------------------------')

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
