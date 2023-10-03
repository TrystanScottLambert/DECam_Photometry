"""
Making postage stamp plots for the paper.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame
from astropy.io import fits
from astropy.visualization import ZScaleInterval
import scipy.ndimage as ndimage

from plot_cumulative import cross_match_to_sexcat
from sex_catalog import SExtractorCat
from plot_postage_stamps import cut_postage_stamp
from zero_points import zero_points
from plotting import end_plot


# Our candidates
INFILE_US = 'candidates_e.txt'
I_CATALOG = '../correct_stacks/N964/i.cat'
Z_CATALOG = '../correct_stacks/N964/z.cat'
N_CATALOG = '../correct_stacks/N964/n964.cat'
I_FITS = '../correct_stacks/N964/i.fits'
Z_FITS = '../correct_stacks/N964/z.fits'
N_FITS = '../correct_stacks/N964/n964.fits'
i_fits = fits.open(I_FITS)
z_fits = fits.open(Z_FITS)
n_fits = fits.open(N_FITS)

SIGMA_I_2_string = '>26.64'
SIGMA_Z_2_string = '>26.58'
SIGMA_I_2 = 26.64
SIGMA_Z_2 = 26.58


os.system('rm postage_stamps/*.png')
def fancy_round(mag) -> str:
    """
    Rounds a given mag value which can be either a flot or a string.
    If it is a string then no rounding is done.
    """
    try:
        val = round(mag, 2)
    except TypeError:
        val = mag
    return str(val)

def replace_non_detections(catalog: DataFrame, replacement_value: str) -> DataFrame:
    """
    Replaces the non detect corrected magnitudes with the replacement value.
    Ususally this replacement value should be the 2 or 3 sigma limiting 
    magnitude. catalog must have the MAG_CORR value calculated.
    """
    catalog.astype({'MAG_CORR': 'str'})
    catalog['SNR'] = (2.5/np.log(10))/catalog['MAGERR_APER']
    catalog.loc[catalog['SNR'] < 2, 'MAG_CORR'] = replacement_value
    return catalog

ra, dec = np.loadtxt(INFILE_US, unpack=True)

i_cat = cross_match_to_sexcat(ra, dec, SExtractorCat(I_CATALOG))
z_cat = cross_match_to_sexcat(ra, dec, SExtractorCat(Z_CATALOG))
n_cat = cross_match_to_sexcat(ra, dec, SExtractorCat(N_CATALOG))

i_cat['MAG_CORR'] = i_cat['MAG_APER'] + zero_points.i_band.mag_correct(1)
z_cat['MAG_CORR'] = z_cat['MAG_APER'] + zero_points.z_band.mag_correct(1)
n_cat['MAG_CORR'] = n_cat['MAG_APER'] + zero_points.n964_band.mag_correct(1)

i_cat['MAG_SNR_KRON'] = (2.5/np.log(10))/i_cat['MAGERR_AUTO']
z_cat['MAG_SNR_KRON'] = (2.5/np.log(10))/z_cat['MAGERR_AUTO']
n_cat['MAG_SNR_KRON'] = (2.5/np.log(10))/n_cat['MAGERR_AUTO']

i_cat = replace_non_detections(i_cat, SIGMA_I_2_string)
z_cat = replace_non_detections(z_cat, SIGMA_Z_2_string)
n_cat = replace_non_detections(n_cat, 'NAS')

i_snr, i_mag = np.array(i_cat['SNR']), np.array(i_cat['MAG_CORR'])
z_snr, z_mag = np.array(z_cat['SNR']), np.array(z_cat['MAG_CORR'])
n_snr, n_mag = np.array(n_cat['SNR']), np.array(n_cat['MAG_CORR'])

for i, ra_val in enumerate(ra):
    z_data = cut_postage_stamp(ra_val, dec[i], z_fits)
    i_data = cut_postage_stamp(ra_val, dec[i], i_fits)
    n_data = cut_postage_stamp(ra_val, dec[i], n_fits)
    z_scale = ZScaleInterval()

    z_data = ndimage.gaussian_filter(z_data, sigma=1.5, order=0)
    i_data = ndimage.gaussian_filter(i_data, sigma=1.5, order=0)
    n_data = ndimage.gaussian_filter(n_data, sigma=1.5, order=0)

    z_min, z_max = z_scale.get_limits(z_data)
    i_min, i_max = z_scale.get_limits(i_data)
    n_min, n_max = z_scale.get_limits(n_data)

    fig = plt.figure()
    ax_i = fig.add_subplot(131)
    ax_n = fig.add_subplot(132)
    ax_z = fig.add_subplot(133)
    ax_i.get_xaxis().set_visible(False)
    ax_i.get_yaxis().set_visible(False)
    ax_z.get_xaxis().set_visible(False)
    ax_z.get_yaxis().set_visible(False)
    ax_n.get_xaxis().set_visible(False)
    ax_n.get_yaxis().set_visible(False)

    ax_i.imshow(i_data, vmin=i_min, vmax=i_max, cmap='gray_r')
    ax_i.scatter(20.5, 20, marker='o', s=800, lw=1, facecolors='none', edgecolors='r')
    ax_n.imshow(n_data, vmin=n_min, vmax=n_max, cmap='gray_r')
    ax_n.scatter(20.5, 20, marker='o', s=800, lw=1, facecolors='none', edgecolors='r')
    ax_z.imshow(z_data, vmin=z_min, vmax=z_max, cmap='gray_r')
    ax_z.scatter(20.5, 20, marker='o', s=800, lw=1, facecolors='none', edgecolors='r')

    ax_i.text(0.75, 0.1, f'{round(i_snr[i], 1)}',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax_i.transAxes, fontsize=20, color='w')
    
    ax_i.text(0.2, 0.1, f'{round(list(i_cat["MAG_SNR_KRON"])[i], 2)}',
              horizontalalignment='center',
              verticalalignment='center',
              transform=ax_i.transAxes, fontsize=20, color='w')

    ax_i.text(0.4, 0.85, fancy_round(i_mag[i]),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax_i.transAxes, fontsize=20, color='w')

    ax_z.text(0.75, 0.1, f'{round(z_snr[i], 1)}',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax_z.transAxes, fontsize=20, color='w')

    ax_z.text(0.4, 0.85, fancy_round(z_mag[i]),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax_z.transAxes, fontsize=20, color='w')

    ax_n.text(0.75, 0.1, f'{round(n_snr[i], 1)}',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax_n.transAxes, fontsize=20, color='w')

    ax_n.text(0.4, 0.85, fancy_round(n_mag[i]),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax_n.transAxes, fontsize=20, color='w')

    end_plot(f'postage_stamps/candidate_{i+1}.png')
    plt.close()
