"""Determining the depths of the of the DECam images."""

import pylab as plt
import numpy as np
from scipy.stats import gaussian_kde
import plotting
from sex_catalog import SExtractorCat
from zero_points import zero_points
from zero_points_cdfs import zero_points_cdfs


HU_VALUES = {
    'i': 27.3,
    'z': 27.0,
    'n964': 25
}

def read_in_sextractor_mags(sextractor_file_name: str):
    """
    Reads in the sextractor mags and uncertainties for the depth plot.
    """
    decam_catalog = SExtractorCat(sextractor_file_name)
    decam_mag = decam_catalog.catalog['MAG_APER'].values
    mag_err = decam_catalog.catalog['MAGERR_APER'].values
    return decam_mag, mag_err


def plot_depth(sextractor_cat_name: str, zpt: float, broadband: str) -> None:
    """
    Makes a mag vs mag_err plot from all the available detections.
    """
    if broadband not in 'izn964':
        raise ValueError('broad band must be either "i", "z", or "n964".')

    mag, mag_err = read_in_sextractor_mags(sextractor_cat_name)
    mag += zpt
    good_values = np.where(mag < 50)[0]
    mag, mag_err = mag[good_values], mag_err[good_values]

    signal_to_noise = (2.5/np.log(10))/mag_err
    signal_to_noise_cut = np.where(signal_to_noise < 30)[0]
    mag, signal_to_noise  = mag[signal_to_noise_cut], signal_to_noise[signal_to_noise_cut]

    # Do some shadding instead of a scatter plot.
    #nbins = 100
    #kernal = gaussian_kde((mag, signal_to_noise))
    #x_i, y_i = np.mgrid[
    #    mag.min():mag.max():nbins*1j, signal_to_noise.min():signal_to_noise.max():nbins*1j]
    #z_i = kernal(np.vstack([x_i.flatten(), y_i.flatten()]))
    plotting.start_plot(x_label='Measured Mags', y_label='SNR')
    #plt.pcolormesh(x_i, y_i, z_i.reshape(x_i.shape), shading='auto')
    plt.scatter(mag, signal_to_noise, s= 0.01, color='k', alpha=0.3)
    plt.axhline(5, ls='--', lw=1.5, color='r', alpha=0.4)
    plt.axvline(HU_VALUES[broadband], ls = ':', color='k')
    plt.xlim(22, 28)
    plotting.end_plot(f'plots/depth_{broadband}.png')
    plt.show()


if __name__ == '__main__':
    INFILE_SEX_I = '../correct_stacks/N964/i_depth.cat'
    INFILE_SEX_Z = '../correct_stacks/N964/z_depth.cat'
    INFILE_SEX_N = '../correct_stacks/N964/n964.cat'
    INFILE_SEX_N_135 = '../correct_stacks/N964/n964_135.cat'
    I_ZPT = zero_points.i_band.mag_correct(1)
    Z_ZPT = zero_points.z_band.mag_correct(1)
    N_ZPT = zero_points.n964_band.mag_correct(1)
    N_ZPT_135 = zero_points.n964_band.mag_correct(1.35/2)


    INFILE_SEX_I_CDFS = '../CDFS_LAGER/i_cdfs_depth.cat'
    INFILE_SEX_Z_CDFS = '../CDFS_LAGER/z_cdfs_depth.cat'
    INFILE_SEX_N_CDFS = '../CDFS_LAGER/n964_cdfs.cat'
    INFILE_SEX_N135_CDFS = '../CDFS_LAGER/n964_135_cdfs.cat'
    I_ZPT_CDFS = zero_points_cdfs.i_band.mag_correct(1)
    Z_ZPT_CDFS = zero_points_cdfs.z_band.mag_correct(1)
    N_ZPT_CDFS = zero_points_cdfs.n964_band.mag_correct(1)
    N_ZPT_135_CDFS = zero_points_cdfs.n964_band.mag_correct(1.35/2)

    plot_depth(INFILE_SEX_I, I_ZPT, 'i')
    plot_depth(INFILE_SEX_Z, Z_ZPT, 'z')
    plot_depth(INFILE_SEX_N, N_ZPT, 'n964')
    #plot_depth(INFILE_SEX_N_135, N_ZPT_135, 'n964')

    #plot_depth(INFILE_SEX_I_CDFS, I_ZPT_CDFS, 'i')
    #plot_depth(INFILE_SEX_Z_CDFS, Z_ZPT_CDFS, 'z')
    #plot_depth(INFILE_SEX_N_CDFS, N_ZPT_CDFS, 'n964')
    #plot_depth(INFILE_SEX_N135_CDFS, N_ZPT_135_CDFS, 'n964')
