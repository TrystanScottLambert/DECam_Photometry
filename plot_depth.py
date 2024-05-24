"""Determining the depths of the of the DECam images."""

import pylab as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
import plotting
from sex_catalog import SExtractorCat
from zero_points import zero_points
from zero_points_cdfs import zero_points_cdfs


#exponential fitting
def straight_line(x, m, b):
    return m*(x) + b


# Define the exponential function to fit, including the (x - k) term
def exp_func(x, a, b, c, k):
    return a * np.exp(b * (x - k)) + c

def exponential_fit(x_data, y_data, p0, x_min=22, x_max=28, k=25):
    # Filter data to include only x values between x_min and x_max
    mask = ~np.isnan(x_data) & ~np.isnan(y_data) & ~np.isinf(x_data) & ~np.isinf(y_data)
    x_data_clean = x_data[mask]
    y_data_clean = y_data[mask]

    # Further filter data to include only x values between x_min and x_max
    range_mask = (x_data_clean >= x_min) & (x_data_clean <= x_max)
    x_data_filtered = x_data_clean[range_mask]
    y_data_filtered = y_data_clean[range_mask]

    # Fit the exponential model to the filtered data
    popt, pcov = curve_fit(exp_func, x_data_filtered, y_data_filtered, p0=p0, maxfev=5000)

    return popt, pcov, x_data_filtered, y_data_filtered


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

def read_data(sextractor_cat_name: str, zpt: float):
    mag, mag_err = read_in_sextractor_mags(sextractor_cat_name)
    mag += zpt
    good_values = np.where(mag < 50)[0]
    mag, mag_err = mag[good_values], mag_err[good_values]
    return mag, mag_err

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

    # exponential fitting
    popt, pcov, x, y = exponential_fit(mag, mag_err, [25, 0.5, 25, 0])
    plt.scatter(x, y)
    x_shit = np.linspace(np.min(x), np.max(x))
    plt.plot(x_shit, exp_func(x_shit, *popt), lw=3, color='r')
    plt.show()
    ##


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
    x_plotting = np.linspace(np.min(mag), np.max(mag), 1000)
    plt.plot(x_plotting, exp_func(x_plotting, *popt))
    plt.axhline(5, ls='--', lw=1.5, color='r', alpha=0.4)
    plt.axvline(HU_VALUES[broadband], ls = ':', color='k')
    plt.xlim(22, 28)
    plotting.end_plot(f'plots/depth_{broadband}.png')
    plt.show()
    return popt, pcov


if __name__ == '__main__':
    INFILE_SEX_I = '../correct_stacks/N964/i_depth.cat'
    INFILE_SEX_Z = '../correct_stacks/N964/z_depth.cat'
    INFILE_SEX_N = '../correct_stacks/N964/n964.cat'
    INFILE_SEX_N_135 = '../correct_stacks/N964/n964_135.cat'
    I_ZPT = zero_points.i_band.mag_correct(1)
    Z_ZPT = zero_points.z_band.mag_correct(1)
    N_ZPT = zero_points.n964_band.mag_correct(1)
    N_ZPT_135 = zero_points.n964_band.mag_correct(1.35/2)

    #INFILE_SEX_I_CDFS = '../CDFS_LAGER/i_cdfs_depth.cat'
    #INFILE_SEX_Z_CDFS = '../CDFS_LAGER/z_cdfs_depth.cat'
    #INFILE_SEX_N_CDFS = '../CDFS_LAGER/n964_cdfs.cat'
    #INFILE_SEX_N135_CDFS = '../CDFS_LAGER/n964_135_cdfs.cat'
    #I_ZPT_CDFS = zero_points_cdfs.i_band.mag_correct(1)
    #Z_ZPT_CDFS = zero_points_cdfs.z_band.mag_correct(1)
    #N_ZPT_CDFS = zero_points_cdfs.n964_band.mag_correct(1)
    #N_ZPT_135_CDFS = zero_points_cdfs.n964_band.mag_correct(1.35/2)

    _, _ = plot_depth(INFILE_SEX_I, I_ZPT, 'i')
    z_popt, _ = plot_depth(INFILE_SEX_Z, Z_ZPT, 'z')
    n_popt, _ = plot_depth(INFILE_SEX_N, N_ZPT, 'n964')
    #plot_depth(INFILE_SEX_N_135, N_ZPT_135, 'n964')

    #plot_depth(INFILE_SEX_I_CDFS, I_ZPT_CDFS, 'i')
    #plot_depth(INFILE_SEX_Z_CDFS, Z_ZPT_CDFS, 'z')
    #plot_depth(INFILE_SEX_N_CDFS, N_ZPT_CDFS, 'n964')
    #plot_depth(INFILE_SEX_N135_CDFS, N_ZPT_135_CDFS, 'n964')
