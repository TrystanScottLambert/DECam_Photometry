"""Main script to check photometry of DECAM images."""

from typing import Tuple
import pylab as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sex_catalog import SExtractorCat

def remove_outliers(array, sigma):
    """returns a mask of values within the given sigma."""
    median, std = np.median(array), np.std(array)
    good_idx = np.where((array >= median-sigma*std) & (array <= median + sigma*std))[0]
    return good_idx

def straight_line(x_array, constant):
    """straight line to fit to data."""
    return constant * np.ones(len(x_array))

def read_in_wide_band(sextractor_file_name: str, panstars_cat_name: str):
    """Reading in the catalogs."""
    decam_catalog = SExtractorCat(sextractor_file_name)
    decam_catalog, pan_cat = decam_catalog.cross_match_with_panstars(panstars_cat_name)
    return decam_catalog, pan_cat

def make_wide_band_mags(
        decam_df: pd.DataFrame, panstars_cat: pd.DataFrame, band: str
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Gets the magnitudes and errors for a wide-band filter."""

    dec_mags = decam_df['MAG_BEST'].values
    dec_mags_uncertainty = decam_df['MAGERR_BEST'].values
    pan_mags = panstars_cat[f'{band}MeanPSFMag'].values
    pan_mags_uncertainty = panstars_cat[f'{band}MeanPSFMagErr'].values

    converted_mags, converted_mags_errors = CONVERT[band](panstars_cat)

    mag_diff, mag_diff_err = calculate_difference_between_mags(
        dec_mags, dec_mags_uncertainty, converted_mags, converted_mags_errors)


    return dec_mags, dec_mags_uncertainty, pan_mags, pan_mags_uncertainty, \
        converted_mags, converted_mags_errors, mag_diff, mag_diff_err

def convert_panstars_i_dec_mags(panstars_cat: pd.DataFrame):
    """Converts the panstars magnitudes into decam magnitudes.
    see: https://des.ncsa.illinois.edu/releases/dr2/dr2-docs/dr2-transformations"""
    r_panstars_mags = panstars_cat['rMeanPSFMag'].values
    i_panstars_mags = panstars_cat['iMeanPSFMag'].values
    i_uncertainties = panstars_cat['iMeanPSFMagErr'].values
    r_uncertainties = panstars_cat['rMeanPSFMagErr'].values

    i_decam = i_panstars_mags - 0.155*(r_panstars_mags - i_panstars_mags) + 0.015

    converted_i_uncertainties = np.hypot(
        i_uncertainties, 0.155*np.hypot(i_uncertainties, r_uncertainties))

    return i_decam, converted_i_uncertainties

def convert_panstars_z_dec_mags(panstars_cat: pd.DataFrame):
    """Converts the panstars magnitudes into decam magnitudes.
    see: https://des.ncsa.illinois.edu/releases/dr2/dr2-docs/dr2-transformations"""
    r_panstars_mags = panstars_cat['rMeanPSFMag'].values
    i_panstars_mags = panstars_cat['iMeanPSFMag'].values
    z_panstars_mags = panstars_cat['zMeanPSFMag'].values
    i_uncertainties = panstars_cat['iMeanPSFMagErr'].values
    r_uncertainties = panstars_cat['rMeanPSFMagErr'].values
    z_uncertainties = panstars_cat['zMeanPSFMagErr'].values

    z_decam = z_panstars_mags - 0.114*(r_panstars_mags - i_panstars_mags) - 0.010

    converted_z_uncertainties = np.hypot(
        z_uncertainties, 0.114*np.hypot(i_uncertainties, r_uncertainties))

    return z_decam, converted_z_uncertainties

CONVERT = {
    'i': convert_panstars_i_dec_mags,
    'z': convert_panstars_z_dec_mags,
}

def calculate_difference_between_mags(
    decam_mags, decam_mags_uncertainties, converted_decam_mags, converted_decam_mags_uncertainties):
    """Working out the y_axis of the plot."""
    difference = decam_mags - converted_decam_mags
    difference_errors = np.hypot(decam_mags_uncertainties, converted_decam_mags_uncertainties)
    return difference, difference_errors

def prepare_plotting_data(sextractor_file_name, panstars_file_name, band):
    """Gets all the necessary variables for plotting ready and removes outliers."""
    decam_catalog, panstars_catalog = read_in_wide_band(sextractor_file_name, panstars_file_name)
    mags = make_wide_band_mags(decam_catalog, panstars_catalog, band = band)
    msk = remove_outliers(mags[-2], sigma=2)
    cleaned_mags = [mag[msk] for mag in mags]
    return cleaned_mags


def plot_direct_comparison(mag_1, mag_1_errors, mag_2, mag_2_errors, x_label = '', y_label = ''):
    """Just plots the two magnitudes against one another"""
    plt.errorbar(mag_1, mag_2, xerr=mag_1_errors, yerr=mag_2_errors, fmt='ko',alpha=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

def plot_difference_fit(mag, diff, mag_err, diff_err, x_lim = None, y_lim = None):
    """Plots the difference plot. (Main plot). Need a difference worked out already.
    left-bottom corner and right-top corner can be included."""
    if x_lim is not None:
        cut1 = np.where((mag <= x_lim[1]) & (mag >= x_lim[0]))[0]
        cut2 = np.where((diff <= y_lim[1]) & (diff >= y_lim[0]))[0]
        cut = np.intersect1d(cut1, cut2)
        mag_cut = mag[cut]
        #mag_err_cut = mag_err[cut]
        diff_cut = diff[cut]
        diff_err_cut = diff_err[cut]

        a_fit, cov = curve_fit(straight_line, mag_cut, diff_cut, sigma=diff_err_cut, absolute_sigma=True)
    else:
        a_fit, cov = curve_fit(straight_line, mag, diff, sigma=diff_err, absolute_sigma=True)
    uncertainties = np.sqrt(np.diag(cov))
    c = a_fit[0]
    plt.title(f'{c} +- {uncertainties[0]}')
    plt.errorbar(mag, diff, xerr=mag_err, yerr=diff_err, fmt='ko', ms=3, capsize=0, alpha=0.5)
    plt.xlabel('decam_mags')
    plt.ylabel('pan_stars_mags converted into decam Mags - decam mags')

    x = np.linspace(np.sort(mag)[0], np.sort(mag)[-1])
    plt.plot(x, straight_line(x, a_fit[0]), ls='--', color='k', lw=2)
    plt.show()


if __name__ == '__main__':
    INFILE_SEX = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/test.cat'
    INFILE_PAN = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_i.csv'
    mags = prepare_plotting_data(INFILE_SEX, INFILE_PAN, band='i')

    plot_direct_comparison(
        mags[0], mags[1],
        mags[2], mags[3],
        x_label='Decam Mags', y_label='Panstars_Mags')

    plot_direct_comparison(
        mags[0], mags[1],
        mags[4], mags[5],
        x_label='Decam Mags', y_label='Panstars Converted into Decam')

    plot_difference_fit(mags[0], mags[6], mags[1], mags[7], x_lim=(-15.208, -11.497), y_lim=(-31.6, -30.44))

    plt.scatter(mags[0], mags[1])
    plt.show()

    INFILE_SEX = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/z/test.cat'
    INFILE_PAN = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_z.csv'
    mags = prepare_plotting_data(INFILE_SEX, INFILE_PAN, band='z')

    plot_direct_comparison(
        mags[0], mags[1],
        mags[2], mags[3],
        x_label='Decam Mags', y_label='Panstars_Mags')

    plot_direct_comparison(
        mags[0], mags[1],
        mags[4], mags[5],
        x_label='Decam Mags', y_label='Panstars Converted into Decam')

    plot_difference_fit(mags[0], mags[6], mags[1], mags[7], x_lim=(-15.9, -12.2), y_lim=(-30.611, -30.455))

    plt.scatter(mags[0], mags[1])
    plt.show()

    # Mag vs Error:
