"""
Module for fitting exponential/splines to SNR-type data.

This is essential for determining things like depth using 
Sextractor uncertainties.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from sex_catalog import SExtractorCat
from zero_points import zero_points, ZeroPoint
from zero_points_cdfs import zero_points_cdfs
from plotting import start_plot, end_plot


def exponential_func(x, a, b):
    """Exponential function for fitting."""
    return a * np.exp(b * x)

def calculate_snr(mag_err: float) -> float:
    """Converts the magnitude error into a snr value."""
    return (2.5/np.log(10))/mag_err

def load_catalog(catalog_name: str, zero_point: ZeroPoint) -> tuple[np.ndarray, np.ndarray]:
    """Loads the catalog for plotting"""
    cat_all = SExtractorCat(catalog_name).catalog
    cat = cat_all[cat_all['MAG_APER'] != 99]
    mag, err = cat['MAG_APER'].values, cat['MAGERR_APER'].values
    cut = np.where((err<50))[0]
    mag, err = mag[cut], err[cut]
    mag_sort = np.argsort(mag)
    mag = mag[mag_sort]
    err = err[mag_sort]
    snr = calculate_snr(err)
    invalid_mask = np.logical_or(np.isnan(err), np.isinf(err))
    valid_indices = np.where(~invalid_mask)
    mag = mag[valid_indices]
    err = err[valid_indices]
    snr = snr[valid_indices]

    cut = np.where(snr<30)[0]
    mag, err, snr = mag[cut], err[cut], snr[cut]
    return mag + zero_point.mag_correct(1), err, snr


if __name__ == '__main__':
    #INFILE = '../correct_stacks/N964/i_depth.cat'
    #zero_point = zero_points.i_band
    #OUTFILE = 'plots/i_band_depth.png'

    INFILE = '../correct_stacks/N964/z_depth.cat'
    zero_point = zero_points.z_band
    OUTFILE = 'plots/z_band_depth.png'
    
    #INFILE = '../correct_stacks/N964/n964.cat'
    #zero_point = zero_points.n964_band
    #OUTFILE = 'plots/n964_band_depth.png'
    
    #INFILE = '../CDFS_LAGER/i_cdfs_depth.cat'
    #zero_point = zero_points_cdfs.i_band
    #OUTFILE = 'plots/i_band_depth_cdfs.png'

    #INFILE = '../CDFS_LAGER/z_cdfs_depth.cat'
    #zero_point = zero_points_cdfs.z_band
    #OUTFILE = 'plots/z_band_depth_cdfs.png'

    #INFILE = '../CDFS_LAGER/n964_cdfs.cat'
    #zero_point = zero_points_cdfs.n964_band
    #OUTFILE = 'plots/n964_band_depth_cdfs.png'


    mag, err, snr = load_catalog(INFILE, zero_point)
    params, cov_matrix = curve_fit(exponential_func, xdata=mag, ydata=snr, maxfev=5000)

    def magnitude_at_snr(snr: float) -> float:
        """Works out the magnitude that results in the given snr"""
        return np.log(snr/params[0])/params[1]


    # Extract the fitting parameters
    a_fit, b_fit = params
    a_err, b_err = np.sqrt(np.diag(cov_matrix))

    # Generating finer x values for smoother plot
    x_finer = np.linspace(mag[0], mag[-1], 1000)

    # Evaluate the fitted exponential curve at the finer x values
    y_fitted = exponential_func(x_finer, a_fit, b_fit)
    y_upper = exponential_func(x_finer, a_fit + 5 *a_err, b_fit + 5 * b_err)
    y_lower = exponential_func(x_finer, a_fit - 5* a_err, b_fit - 5 * b_err)


    # Plotting the original data and the fitted exponential curve

    start_plot('Magnitude', 'SNR')
    plt.scatter(mag, snr, s=1)
    plt.plot(x_finer, y_fitted, label='Exponential Fit', color='red')
    plt.fill_between(x_finer, y_lower, y_upper, color='red', alpha=0.2, label='Uncertainty')
    plt.xlim(left=np.min(mag))
    plt.ylim(bottom=0)
    end_plot(OUTFILE)

    #printing out the depths
    for snr in range(1, 6):
        print(f'for SNR = {snr} depth is: {magnitude_at_snr(snr)}')
