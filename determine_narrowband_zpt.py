"""
Determining the narrow band zpt.
Have to use a different method than the broad bands because
DECAM doesn't have a calibration method. Using the method suggested
to us via private communication from ...
"""

import pylab as plt
import numpy as np
from scipy import odr
from sex_catalog import SExtractorCat
import plotting
from k_constant import calculate_k_constant_mag

def straight_line(parameters, x_val):
    """Straight line for fitting to data"""
    alpha, zpt = parameters
    return x_val * alpha + zpt

# bottom and upper limits are there to exclude outliers. The mc values were determined through trial and error.
def bottom_limit(x_val):
    """Lower limit to remove outlier values when making the fit."""
    m = 1
    c = 28.8
    return m * x_val + c

def upper_limit(x_val):
    """Upper limit to remove the upper outliers."""
    m = 1
    c = 29.1
    return m * x_val + c

def is_point_in_limits(x_val:float, y_val:float):
    """Determines if a cartesian point is withint the bottom_limit and upper_limit"""
    upper_lim = upper_limit(x_val)
    bottom_lim = bottom_limit(x_val)
    if y_val < upper_lim and y_val > bottom_lim:
        return True
    else:
        return False

def index_in_limits(x_array, y_array):
    """Returns an array of indicies of all the values within the limits"""
    idx = [i for i in range(len(x_array)) if is_point_in_limits(x_array[i], y_array[i])]
    return np.array(idx)

if __name__ == '__main__':
    INFILE_SEX = '../correct_stacks/N964/n964.cat'
    INFILE_Z = '../PANSTARS/PANSTARS_z.csv'
    INFILE_Y = '../PANSTARS/PANSTARS_y.csv'
    NARROW_BAND_SEEING = 1.32 # arcseconds
    NARROW_BAND_APERTURES = 1. # arcseconds

    decam_all = SExtractorCat(INFILE_SEX)
    decam_catalog_z, pan_cat_z = decam_all.cross_match_with_panstars(INFILE_Z)
    decam_catalog_y, pan_cat_y = decam_all.cross_match_with_panstars(INFILE_Y)

    z = pan_cat_z['zMeanPSFMag'].values
    y = pan_cat_z['yMeanPSFMag'].values
    z_err = pan_cat_z['zMeanPSFMagErr'].values
    y_err = pan_cat_z['zMeanPSFMagErr'].values
    n964 = decam_catalog_z['MAG_AUTO'].values
    n964_err = decam_catalog_z['MAGERR_AUTO'].values
    cut = np.where(y != -999)[0]
    z, y, n964, z_err, y_err, n964_err = z[cut], y[cut], n964[cut], z_err[cut], y_err[cut], n964_err[cut]

    delta_pan = z - y
    delta_decam = z - n964
    delta_pan_err = np.hypot(z_err, y_err)
    delta_decam_err = np.hypot(z_err, n964_err)
    cut_x = np.where((delta_pan > -0.2) & (delta_pan < 0.6))[0]
    cut_y = np.where((delta_decam > 28.75) & (delta_decam < 29.752))
    final_cut = np.intersect1d(cut_x, cut_y)
    fit_cut = index_in_limits(delta_pan[final_cut], delta_decam[final_cut])

    straight_line_model = odr.Model(straight_line)
    data = odr.RealData(delta_pan[final_cut][fit_cut], delta_decam[final_cut][fit_cut],
                        sx=delta_pan[final_cut][fit_cut], sy=delta_decam[final_cut][fit_cut])

    _odr = odr.ODR(data, straight_line_model, beta0=[0.81, 28.9])
    out = _odr.run()
    popt = out.beta
    perr = out.sd_beta
    print("fit parameter 1-sigma error")
    print("———————————–")
    for i, pop in enumerate(popt):
        print(f'{pop} +- {perr[i]}')

    k_const = calculate_k_constant_mag(NARROW_BAND_APERTURES, NARROW_BAND_SEEING)
    zpt_prime = popt[1] - k_const
    print('The zpt prime is: ', zpt_prime)

    NSTD = 5.
    popt_up = popt + NSTD * perr
    popt_dw = popt - NSTD * perr
    x_fit = np.linspace(-0.3, 0.7, 100)
    fit = straight_line(popt, x_fit)
    fit_up = straight_line(popt_up, x_fit)
    fit_dw= straight_line(popt_dw, x_fit)

    fig = plotting.start_plot(x_label = 'z - y', y_label = 'N964 - z')
    plt.errorbar(delta_pan[final_cut], delta_decam[final_cut], xerr=delta_pan_err[final_cut],
                  yerr = delta_decam_err[final_cut], fmt = 'ko', alpha=0.3, capsize=0, zorder=0)

    # To check that the limits are in the right place
    #plt.plot(x_fit, bottom_limit(x_fit), color='b')
    #plt.plot(x_fit, upper_limit(x_fit), color='b')
    plt.plot(x_fit, fit,'r', lw=2)
    plt.fill_between(x_fit, fit_up, fit_dw, alpha=.5, color='r')
    plt.xlim(-0.3, 0.7)
    plotting.end_plot('plots/N964_zpt_calculation.png')

    N964_DIR = '../correct_stacks/N964/'
    N964_SEX_CAT = N964_DIR + 'n964.cat'
    N964_DECAM_CAT = N964_DIR + 'c4d_210831_050404_osj_N964_vik2.fits.fz'

    n964_sex_cat = SExtractorCat(N964_SEX_CAT)
    n964_sex_cat.catalog = n964_sex_cat.catalog[n964_sex_cat.catalog['MAG_AUTO'] < 3]
    n964_sex_cat.to_region_file(N964_DECAM_CAT, N964_DIR + 'N964.reg')
