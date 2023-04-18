"""
Calculate the zero points of the i-band and z-band
"""

from typing import Tuple
import numpy as np
import pylab as plt
from scipy.optimize import curve_fit
import pandas as pd
from sex_catalog import SExtractorCat


def remove_outliers(array, sigma):
    """returns a mask of values within the given sigma."""
    median, std = np.median(array), np.std(array)
    good_idx = np.where((array >= median-sigma*std) & (array <= median + sigma*std))[0]
    return good_idx

def read_in_wide_band(sextractor_file_name: str, panstars_cat_name: str):
    """Reading in the catalogs."""
    decam_catalog = SExtractorCat(sextractor_file_name)
    decam_catalog, pan_cat = decam_catalog.cross_match_with_panstars(panstars_cat_name)
    return decam_catalog, pan_cat

def get_broad_band_mags(decam_df: pd.DataFrame, panstars_cat: pd.DataFrame, band: str) -> Tuple:
    """Gets the magnitudes and errors for a wide-band filter."""

    dec_mags = decam_df['MAG_AUTO'].values
    dec_mags_uncertainty = decam_df['MAGERR_AUTO'].values
    pan_mags = panstars_cat[f'{band}MeanPSFMag'].values
    pan_mags_uncertainty = panstars_cat[f'{band}MeanPSFMagErr'].values

    converted_mags, converted_mags_errors = CONVERT[band](panstars_cat)

    return dec_mags, dec_mags_uncertainty, converted_mags, converted_mags_errors, pan_mags, pan_mags_uncertainty


def convert_panstars_i_dec_mags(panstars_cat: pd.DataFrame):
    """Converts the panstars magnitudes into decam magnitudes.
    see: https://des.ncsa.illinois.edu/releases/dr2/dr2-docs/dr2-transformations"""
    r_panstars_mags = panstars_cat['rMeanPSFMag'].values
    i_panstars_mags = panstars_cat['iMeanPSFMag'].values
    i_uncertainties = panstars_cat['iMeanPSFMagErr'].values
    r_uncertainties = panstars_cat['rMeanPSFMagErr'].values

    i_decam = i_panstars_mags - 0.155 * (r_panstars_mags - i_panstars_mags) + 0.015
    converted_i_uncertainties = np.hypot(i_uncertainties, 0.155*np.hypot(i_uncertainties, r_uncertainties))

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

    z_decam = z_panstars_mags - 0.114 * (r_panstars_mags - i_panstars_mags) - 0.010
    converted_z_uncertainties = np.hypot(z_uncertainties, 0.114*np.hypot(i_uncertainties, r_uncertainties))

    return z_decam, converted_z_uncertainties


CONVERT = {
    'i': convert_panstars_i_dec_mags,
    'z': convert_panstars_z_dec_mags,
}


def prepare_plotting_data(sextractor_file_name, panstars_file_name, band):
    """Gets all the necessary variables for plotting ready and removes outliers."""
    decam_catalog, panstars_catalog = read_in_wide_band(sextractor_file_name, panstars_file_name)
    mags = get_broad_band_mags(decam_catalog, panstars_catalog, band=band)
    msk = remove_outliers(mags[-2], sigma=2)
    cleaned_mags = [mag[msk] for mag in mags]
    return cleaned_mags


def straight_line(x_array, intercept):
    """y=mx +c with m=1. Line for fitting to the zpt plot."""
    return x_array + intercept

class BroadBand:
    """
    Broad Band class for sorting the different magnitudes
    and determine the zpt.
    """
    def __init__(self, panstars_cat_name: str, sextractor_cat_name: str, broadband: str):
        if broadband not in 'iz':
            raise ValueError('broadband needs to be either "i" or "z".')
        self.sextractor_cat_name = sextractor_cat_name
        mags = prepare_plotting_data(sextractor_cat_name, panstars_cat_name, broadband)
        good_values = np.where(mags[2] < 100)[0] # Obviously not physical if this isn't met.
        self.measured_mags = mags[0][good_values]
        self.measured_mags_err = mags[1][good_values]
        self.expected_mags = mags[2][good_values]
        self.expected_mags_err = mags[3][good_values]

    def plot_zpt(self):
        """Plots the measured mags vs the expected mags."""
        plt.errorbar(
        self.measured_mags, self.expected_mags,
        xerr=self.measured_mags_err, yerr=self.expected_mags_err,
        fmt='ro', alpha=0.3)
        x_fit, fit, fit_up, fit_down = self.fit_straight_line()
        plt.plot(x_fit, fit, ls='--', color='k', lw=3, zorder=1)
        plt.fill_between(x_fit, fit_up, fit_down, alpha=.5, color='r')
        plt.show()

    def fit_straight_line(self) -> Tuple:
        """Fitting y=x+c line to the data to determine c"""
        zpt, zpt_err = self.zero_point
        nstd = 5.
        popt_up = zpt + nstd * zpt_err
        popt_dw = zpt - nstd * zpt_err
        x_fit = np.linspace(np.sort(self.measured_mags)[0], np.sort(self.measured_mags)[-1])
        fit = straight_line(x_fit, zpt)
        fit_up = straight_line(x_fit, popt_up)
        fit_dw= straight_line(x_fit, popt_dw)
        return x_fit, fit, fit_up, fit_dw

    @property
    def zero_point(self) -> Tuple[float, float]:
        a_fit, cov = curve_fit(
            straight_line, self.measured_mags, self.expected_mags,
            sigma=self.expected_mags_err, absolute_sigma=True)
        uncertainties = np.sqrt(np.diag(cov))
        return a_fit[0], uncertainties[0]


if __name__ == '__main__':
    INFILE_SEX_I = '../correct_stacks/N964/i.cat'
    INFILE_PAN_I = '../PANSTARS/PANSTARS_i.csv'
    INFILE_SEX_Z = '../correct_stacks/N964/z.cat'
    INFILE_PAN_Z = '../PANSTARS/PANSTARS_z.csv'

    i_band = BroadBand(INFILE_PAN_I, INFILE_SEX_I, 'i')
    i_band.plot_zpt()

    z_band = BroadBand(INFILE_PAN_Z, INFILE_SEX_Z, 'z')
    z_band.plot_zpt()
