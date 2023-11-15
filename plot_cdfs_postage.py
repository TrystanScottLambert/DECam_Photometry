"""
Making postage stamp plot for the cdfs galaxies.
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
from zero_points_cdfs import zero_points_cdfs
from plotting import end_plot


# Our candidates
#INFILE_US = 'candidates_cdfs_e.txt'
INFILE_US = 'candidates_cdfs_e.txt'
I_CATALOG = '../CDFS_LAGER/i_cdfs.cat'
Z_CATALOG = '../CDFS_LAGER/z_cdfs.cat'
N_CATALOG = '../CDFS_LAGER/n964_cdfs.cat'
I_FITS = '../CDFS_LAGER/CDFS_i.fits'
Z_FITS = '../CDFS_LAGER/CDFS_z.fits'
N_FITS = '../CDFS_LAGER/CDFS_NB.fits'

i_fits = fits.open(I_FITS)
z_fits = fits.open(Z_FITS)
n_fits = fits.open(N_FITS)

SIGMA_I_2_string = '>28.10'
SIGMA_Z_2_string = '>27.73'
N964_INFINITY = -999
SIGMA_I_2 = 28.10
SIGMA_Z_2 = 27.73

if INFILE_US == 'candidates_true_cdfs.txt':
    os.system('rm postage_stamps_cdfs/true/*.png')
else:
    os.system('rm postage_stamps_cdfs/degraded/*.png')

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

def replace_non_detections(catalog: DataFrame, replacement_value: str, limit: float) -> DataFrame:
    """
    Replaces the non detect corrected magnitudes with the replacement value.
    Ususally this replacement value should be the 2 or 3 sigma limiting 
    magnitude. catalog must have the MAG_CORR value calculated.
    """
    catalog.astype({'MAG_CORR': 'str'})
    catalog['SNR'] = (2.5/np.log(10))/catalog['MAGERR_APER']
    catalog.loc[catalog['SNR'] < limit, 'MAG_CORR'] = replacement_value
    return catalog

ra, dec = np.loadtxt(INFILE_US, unpack=True)

i_cat = cross_match_to_sexcat(ra, dec, SExtractorCat(I_CATALOG))
z_cat = cross_match_to_sexcat(ra, dec, SExtractorCat(Z_CATALOG))
n_cat = cross_match_to_sexcat(ra, dec, SExtractorCat(N_CATALOG))

i_cat['MAG_CORR'] = i_cat['MAG_APER'] + zero_points_cdfs.i_band.mag_correct(1)
z_cat['MAG_CORR'] = z_cat['MAG_APER'] + zero_points_cdfs.z_band.mag_correct(1)
n_cat['MAG_CORR'] = n_cat['MAG_APER'] + zero_points_cdfs.n964_band.mag_correct(1)

i_cat['MAG_SNR_KRON'] = (2.5/np.log(10))/i_cat['MAGERR_AUTO']
z_cat['MAG_SNR_KRON'] = (2.5/np.log(10))/z_cat['MAGERR_AUTO']
n_cat['MAG_SNR_KRON'] = (2.5/np.log(10))/n_cat['MAGERR_AUTO']

i_cat = replace_non_detections(i_cat, SIGMA_I_2_string, 2)
z_cat = replace_non_detections(z_cat, SIGMA_Z_2_string, 2)
n_cat = replace_non_detections(n_cat, 'NAS', 2)

i_snr, i_mag = np.array(i_cat['SNR']), np.array(i_cat['MAG_CORR'])
z_snr, z_mag = np.array(z_cat['SNR']), np.array(z_cat['MAG_CORR'])
n_snr, n_mag = np.array(n_cat['SNR']), np.array(n_cat['MAG_CORR'])

for i, ra_val in enumerate(ra):
    z_data = cut_postage_stamp(ra_val, dec[i], z_fits, hdu_number=0)
    i_data = cut_postage_stamp(ra_val, dec[i], i_fits, hdu_number=0)
    n_data = cut_postage_stamp(ra_val, dec[i], n_fits, hdu_number=0)

    z_scale = ZScaleInterval()
    z_min, z_max = z_scale.get_limits(z_data)
    i_min, i_max = z_scale.get_limits(i_data)
    n_min, n_max = z_scale.get_limits(n_data)

    #z_data = ndimage.gaussian_filter(z_data, sigma=(1, 1), order=0)
    #i_data = ndimage.gaussian_filter(i_data, sigma=(1, 1), order=0)
    #n_data = ndimage.gaussian_filter(n_data, sigma=(1, 1), order=0)

    fig = plt.figure()
    ax_i = fig.add_subplot(131)
    ax_n = fig.add_subplot(132)
    ax_z = fig.add_subplot(133)
    ax_i.set_xticks([])
    ax_i.set_yticks([])
    ax_z.set_xticks([])
    ax_z.set_yticks([])
    ax_n.set_xticks([])
    ax_n.set_yticks([])

    ax_i.imshow(i_data, vmin=i_min, vmax=i_max, cmap='gray_r')
    #ax_i.scatter(20.5, 20, marker='o', s=800, lw=1, facecolors='none', edgecolors='r')
    ax_i.plot([10,15],[20, 20], lw=2, color='r')
    ax_i.plot([40-10,40-15],[20, 20], lw=2, color='r')
    ax_i.plot([20,20],[10, 15], lw=2, color='r')
    ax_i.plot([20,20],[40-10, 40-15], lw=2, color='r')

    ax_n.imshow(n_data, vmin=n_min, vmax=n_max, cmap='gray_r')
    #ax_n.scatter(20.5, 20, marker='o', s=800, lw=1, facecolors='none', edgecolors='r')
    ax_n.plot([10,15],[20, 20], lw=2, color='r')
    ax_n.plot([40-10,40-15],[20, 20], lw=2, color='r')
    ax_n.plot([20,20],[10, 15], lw=2, color='r')
    ax_n.plot([20,20],[40-10, 40-15], lw=2, color='r')

    ax_z.imshow(z_data, vmin=z_min, vmax=z_max, cmap='gray_r')
    #ax_z.scatter(20.5, 20, marker='o', s=800, lw=1, facecolors='none', edgecolors='r')
    ax_z.plot([10,15],[20, 20], lw=2, color='r')
    ax_z.plot([40-10,40-15],[20, 20], lw=2, color='r')
    ax_z.plot([20,20],[10, 15], lw=2, color='r')
    ax_z.plot([20,20],[40-10, 40-15], lw=2, color='r')


    ax_i.text(0.75, 0.1, f'{round(i_snr[i], 1)}',
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax_i.transAxes, fontsize=20, color='w')

    ax_i.text(0.2, 0.1, f'{round(list(i_cat["MAGERR_AUTO"])[i], 2)}',
              horizontalalignment='center',
              verticalalignment='center',
              transform=ax_i.transAxes, fontsize=20, color='w')

    if i == len(ra) -1:
        ax_i.set_xlabel('i', fontsize=12)
        ax_z.set_xlabel('z', fontsize=12)
        ax_n.set_xlabel('NB964', fontsize=12)

    ax_i.set_ylabel(f'CDFS-{i+1}', fontsize=12)

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

    if INFILE_US == 'candidates_true_cdfs.txt':
        end_plot(f'postage_stamps_cdfs/true/candidate_cdfs_true{i+1}.png')
    else:
        end_plot(f'postage_stamps_cdfs/degraded/candidate_cdfs{i+1}.png')
    plt.close()
