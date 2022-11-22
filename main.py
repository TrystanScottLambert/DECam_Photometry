"""Main script to check photometry of DECAM images."""

import pylab as plt
import numpy as np
from scipy.optimize import curve_fit
#from decam import DecamImage
from sex_catalog import SExtractorCat


def straight_line(x, c, m):
    """straight line to fit to data."""
    return m*x + c

#INFILE = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/c4d_211021_003940_osj_i_vik1_skysubtracted_flux.fits.fz'
INFILE  = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/test.cat'
PAN_CAT = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_i.csv'
#i_band = DecamImage(INFILE)
#decam_cat, pan_cat = i_band.cross_match_with_panstars(PAN_CAT)
#decam_mags = decam_cat['mag'] + 25.128

decam_catalog = SExtractorCat(INFILE)
decam_cat, pan_cat = decam_catalog.cross_match_with_panstars(PAN_CAT)
decam_mags = decam_cat['MAG_BEST']
decam_mags_uncertainty = decam_cat['MAGERR_BEST']
pan_mags = pan_cat['iMeanPSFMag']
pan_mags_uncertainty = pan_cat['iMeanPSFMagErr']

converted_dec_mags = pan_mags - 0.155 * (pan_cat['rMeanPSFMag'] - pan_mags) + 0.015
converted_dec_mags_error = np.hypot(pan_mags_uncertainty, 0.155*np.hypot(pan_mags_uncertainty, pan_cat['rMeanPSFMagErr']))
cut = np.where(converted_dec_mags < 30)[0]

difference = np.array(list(decam_mags)) - np.array(list(converted_dec_mags))
difference_errors = np.hypot(np.array(list(converted_dec_mags_error)), np.array(list(decam_mags_uncertainty)))

plt.errorbar(decam_mags, pan_mags, xerr=decam_mags_uncertainty,
             yerr=pan_mags_uncertainty, fmt='ro', ms=3, capsize=5,ecolor='k')
plt.xlabel('decam_mags')
plt.ylabel('pan_stars_mags')
plt.show()

plt.errorbar(
    np.array(decam_mags)[cut], 
    np.array(converted_dec_mags)[cut], 
    xerr=np.array(decam_mags_uncertainty)[cut],
    yerr=np.array(converted_dec_mags_error)[cut], 
    fmt='ro', ms=3, capsize=5, ecolor='k')

plt.xlabel('decam_mags')
plt.ylabel('pan_stars_mags converted into decam Mags')
plt.show()

x_data = np.array(decam_mags)[cut]
y_data = np.array(difference)[cut]
delta_y = np.array(difference_errors)[cut]
delta_x = np.array(decam_mags_uncertainty)[cut]
cut = np.where((x_data <= -12) & (x_data >= -14))[0]
cut2 = np.where((y_data < -30.512) & (y_data > -31.41))[0]
cut = np.intersect1d(cut, cut2)
a_fit, cov = curve_fit(
    straight_line, x_data[cut], y_data[cut],sigma=delta_y[cut],absolute_sigma=True)
uncertainties = np.sqrt(np.diag(cov))
c = a_fit[0]
m = a_fit[1]
plt.errorbar(x_data, y_data, xerr=delta_x,
             yerr=delta_y, fmt='ro', ms=3, capsize=5, ecolor='k')
plt.xlabel('decam_mags')
plt.ylabel('pan_stars_mags converted into decam Mags - decam mags')

x = np.linspace(np.sort(x_data)[0], np.sort(x_data)[-1])
plt.plot(x, straight_line(x, *a_fit), ls='--', color='k', lw=2)
plt.ylim(-31.41, -30.512)
plt.show()
