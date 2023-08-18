"""
Script to plot the cummulative plot using the selection criteria from
Eduardo and Chira in 2013 and 2017.
"""
import numpy as np
import pylab as plt
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.io import fits

from plot_cumulative import cross_match_to_sexcat
from sex_catalog import SExtractorCat
from zero_points import zero_points
from zero_points_cdfs import zero_points_cdfs
from plotting import start_plot, end_plot


RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24) * u.deg
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1 *u.deg
qso_position = SkyCoord(ra = RA_QSO, dec = DEC_QSO)
cdfs_hdu = fits.open('CDFS_MASK.fits')
decam_hdu = fits.open('DECAm_MASK.fits')
pix_scale = cdfs_hdu[0].header['PC2_2'] * 3600 #arcseconds
pix_decam = decam_hdu[0].header['PC2_2'] * 3600

CDFS_AREA = len(np.where(cdfs_hdu[0].data == 1)[0]) * (pix_scale**2) # Arcseconds
DECAM_AREA = len(np.where(decam_hdu[0].data == 1)[0]) * (pix_decam**2)
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
inner = len(np.where(results == True)[0])
outer = len(np.where(results == False)[0])
inner_area = 10176726.0573 # from plot_radial_distribution
OUTER_AREA = 27036703.3248 # square arcseconds.


#Error propagation
def propagate_ratio(
        numerator: float, denominator: float, u_numerator: float, u_denominator: float) -> float:
    """Propagates the uncertainty based on the numerator and denominator uncertainties"""
    ratio = numerator/denominator
    uncertainty = ratio * np.hypot(u_numerator/numerator, u_denominator/denominator)
    return ratio, uncertainty

cdfs_density =  len(ra_cdfs)/CDFS_AREA
our_density = len(ra_us)/DECAM_AREA
cdfs_uncertainty = np.sqrt(len(ra_cdfs))/CDFS_AREA
our_uncertainty = np.sqrt(len(ra_us))/DECAM_AREA
inner_density = inner/inner_area
outer_density  = outer/OUTER_AREA
inner_uncertainty = np.sqrt(inner)/inner_area
outer_uncertainty = np.sqrt(outer)/OUTER_AREA
inner_ratio = propagate_ratio(inner_density, cdfs_density, inner_uncertainty, cdfs_uncertainty)
outer_ratio = propagate_ratio(outer_density, cdfs_density, outer_uncertainty, cdfs_uncertainty)
ratio, uncertainty = propagate_ratio(our_density, cdfs_density, our_uncertainty, cdfs_uncertainty)
main_result_ratio = propagate_ratio(
    outer_density, inner_density, outer_uncertainty, inner_uncertainty)

#Printing values
print('CDFS density [arcseconds ^{-2}]: ', cdfs_density, '+-', cdfs_uncertainty)
print('our density [arcseconds ^{-2})]: ', our_density, '+-', our_uncertainty)
print('Ratio: ',  ratio, '+-', uncertainty)
print('-----------------------------')
print('Inner 75cMpc is:')
print('Density[arcseconds ^{-2})]: ', inner_density, '+-', inner_uncertainty)
print('Ratio: ', inner_ratio[0], '+-', inner_ratio[1])
print('-----------------------------')
print('Outer 75 cMpc')
print('Density[arcseconds ^{-2})]: ', outer_density, '+-', outer_uncertainty)
print('Ratio: ', outer_ratio[0], '+-', outer_ratio[1])
print('-----------------------------')
print('Ratio of inner density and outer density:')
print('Ratio: ', main_result_ratio[0], '+-', main_result_ratio[1])

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
