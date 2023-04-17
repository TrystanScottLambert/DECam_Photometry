"""
Calculate the zero points of the i-band and z-band
"""


from typing import Tuple
import numpy as np
import pylab as plt
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

def make_wide_band_mags(decam_df: pd.DataFrame, panstars_cat: pd.DataFrame, band: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
    mags = make_wide_band_mags(decam_catalog, panstars_catalog, band=band)
    msk = remove_outliers(mags[-2], sigma=2)
    cleaned_mags = [mag[msk] for mag in mags]
    return cleaned_mags


if __name__ == '__main__':
    INFILE_SEX_I = '../correct_stacks/N964/i.cat'
    INFILE_PAN_I = '../PANSTARS/PANSTARS_i.csv'
    INFILE_SEX_Z = '../correct_stacks/N964/z.cat'
    INFILE_PAN_Z = '../PANSTARS/PANSTARS_z.csv'

    i_mags = prepare_plotting_data(INFILE_SEX_I, INFILE_PAN_I, 'i')
    plt.errorbar(i_mags[0], i_mags[2], xerr=i_mags[1], yerr=i_mags[3], fmt='ro')
    plt.show()

    z_mags = prepare_plotting_data(INFILE_SEX_Z, INFILE_PAN_Z,'z')
