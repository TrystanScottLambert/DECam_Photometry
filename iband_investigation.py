""" 
Investigating why we are finding i-band detections when there shouldn't be any.
Some of the i-band detections have very obvious i-band emission but have very low
SNR so this script is for generating tools to try figure out why this is the case.
"""

import numpy as np 
from pandas import DataFrame
import pylab as plt

from identify_candidates import calculate_snr, write_region_file
from sex_catalog import SExtractorCat

def add_snr(catalog: DataFrame) -> DataFrame:
    """Updates the data frame with the two types of SNR vals"""
    catalog['SNR_KRON'] = calculate_snr(catalog['MAGERR_AUTO'])
    catalog['SNR_APP'] = calculate_snr(catalog['MAGERR_APER'])
    return catalog


def find_both_non_detections(catalog: DataFrame) -> np.ndarray:
    """
    Finds the case where both the SNR KRON and the SNR APP
    is detect low signal to noise sources.
    """
    cut  = np.where((catalog['SNR_KRON'] < 2) & (catalog['SNR_APP'] < 2))[0]
    return cut

def find_only_kron_non_detections(catalog: DataFrame) -> np.ndarray:
    """Finds the cases where the kron snr is less than two 
    but the apperture snr is high"""
    cut = np.where((catalog['SNR_KRON'] < 2) & (catalog['SNR_APP'] > 2))[0]
    return cut

def find_only_aper_non_detections(catalog: DataFrame) -> np.ndarray:
    """
    Finds the cases where the souce is detected in the kron aperture 
    but not in the apperture apperture.
    """
    cut = np.where((catalog['SNR_KRON'] > 2) & (catalog['SNR_APP'] < 2))[0]
    return cut

def get_positions(catalog: DataFrame, indicies: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Works out the ra and dec array for the given indiciees"""
    ra = np.array(catalog['ALPHAPEAK_J2000'])
    dec = np.array(catalog['DELTAPEAK_J2000'])
    return ra[indicies], dec[indicies]


if __name__ == '__main__':
    cdfs_n_cat = add_snr(SExtractorCat('../CDFS_LAGER/n964_cdfs.cat').catalog)
    cdfs_i_cat = add_snr(SExtractorCat('../CDFS_LAGER/i_cdfs.cat').catalog)
    cdfs_z_cat = add_snr(SExtractorCat('../CDFS_LAGER/z_cdfs.cat').catalog)

    us_n_cat = add_snr(SExtractorCat('../correct_stacks/N964/n964.cat').catalog)
    us_i_cat = add_snr(SExtractorCat('../correct_stacks/N964/i.cat').catalog)
    us_z_cat = add_snr(SExtractorCat('../correct_stacks/N964/z.cat').catalog)


    cdfs_both_non_i = find_both_non_detections(cdfs_i_cat)
    cdfs_only_non_app_i = find_only_aper_non_detections(cdfs_i_cat)
    cdfs_only_non_kron_i = find_only_kron_non_detections(cdfs_i_cat)


    write_region_file(*get_positions(cdfs_i_cat, cdfs_both_non_i), 'non_detections_in_both_i_cdfs.reg')
    write_region_file(*get_positions(cdfs_i_cat, cdfs_only_non_app_i), 'non_detections_in_app_i_cdfs.reg')
    write_region_file(*get_positions(cdfs_i_cat, cdfs_only_non_kron_i), 'non_detections_in_kron_i_cdfs.reg')
