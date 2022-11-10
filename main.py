"""Main script to check photometry of DECAM images."""

import pylab as plt
from decam import DecamImage

INFILE = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/c4d_211021_003940_osj_i_vik1_skysubtracted_flux.fits.fz'
PAN_CAT = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_i.csv'
i_band = DecamImage(INFILE)
decam_cat, pan_cat = i_band.cross_match_with_panstars(PAN_CAT)
decam_mags = decam_cat['mag'] + 25.128
pan_mags = pan_cat['iMeanPSFMag']

converted_dec_mags = pan_mags - 0.155 * (pan_cat['rMeanPSFMag'] - pan_mags) + 0.015
plt.scatter(decam_mags, converted_dec_mags)
plt.xlabel('decam_mags')
plt.ylabel('pan_stars_mags')
plt.show()
