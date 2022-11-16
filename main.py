"""Main script to check photometry of DECAM images."""

import pylab as plt
import numpy as np
#from decam import DecamImage
from sex_catalog import SExtractorCat

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

plt.errorbar(decam_mags, pan_mags, xerr=decam_mags_uncertainty,
             yerr=pan_mags_uncertainty, fmt='ro', ms=3, capsize=5)
plt.xlabel('decam_mags')
plt.ylabel('pan_stars_mags')
plt.show()

plt.errorbar(np.array(decam_mags)[cut], np.array(converted_dec_mags)[cut], xerr=np.array(decam_mags_uncertainty)[cut], 
             yerr=np.array(converted_dec_mags_error)[cut], fmt='ro', ms=3, capsize=5)
plt.xlabel('decam_mags')
plt.ylabel('pan_stars_mags converted into decam Mags')
plt.show()