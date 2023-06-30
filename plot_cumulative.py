"""
Making the cummulative plot to show the overdensity.
"""

import numpy as np
import pylab as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from pandas import DataFrame

from sex_catalog import SExtractorCat
from zero_points import zero_points, ZeroPoints
from zero_points_cdfs import zero_points_cdfs
from plotting import start_plot, end_plot

N135_LIM = 25.10
N_LIM = 24.66
Z_LIM = 25.58
depth = -np.arange(1, -0.01, -0.01)
CDFS_AREA = 20385305 # Arcseconds
DECAM_AREA = 2.87e7 # Arcseconds
IMACS_AREA = 5.01046e6 # Arcseconds

def cross_match_to_sexcat(
        ra_array: np.ndarray, dec_array:np.ndarray, sex_catalog: SExtractorCat
        ) -> DataFrame:
    """Reduces the sexcatalog to the cross matched values"""
    sex_ra, sex_dec = np.array(
        sex_catalog.catalog['ALPHAPEAK_J2000']), np.array(sex_catalog.catalog['DELTAPEAK_J2000'])
    catalog = SkyCoord(ra = sex_ra*u.deg, dec = sex_dec*u.deg)
    candidate_coords = SkyCoord(ra = ra_array * u.deg, dec = dec_array*u.deg)
    idx, _, _ = candidate_coords.match_to_catalog_sky(catalog)
    return sex_catalog.catalog.iloc[idx]

def add_ab_mags(n_cat: DataFrame, n135_cat: DataFrame, z_cat: DataFrame, zero_point: ZeroPoints) -> None:
    """Adds the correct AB mags to the catalog."""
    n_cat['MAG_CORR'] = n_cat['MAG_APER'] + zero_point.n964_band.mag_correct(1)
    n135_cat['MAG_CORR'] = n135_cat['MAG_APER'] + zero_point.n964_band.mag_correct(1.35/2)
    z_cat['MAG_CORR'] = z_cat['MAG_APER'] + zero_point.z_band.mag_correct(1)

def prepare_data(n_cat_name: str, n135_cat_name: str, z_cat_name:str, candidates: str, zero_point:ZeroPoints):
    """Only select the candidate catalogs."""
    n_cat = SExtractorCat(n_cat_name)
    n135_cat = SExtractorCat(n135_cat_name)
    z_cat = SExtractorCat(z_cat_name)

    ra, dec = np.loadtxt(candidates, unpack=True)
    
    matched_n_cat = cross_match_to_sexcat(ra, dec, n_cat)
    matched_n135_cat = cross_match_to_sexcat(ra, dec, n135_cat)
    matched_z_cat = cross_match_to_sexcat(ra, dec, z_cat)

    add_ab_mags(matched_n_cat, matched_n135_cat, matched_z_cat, zero_point)
    return matched_n_cat, matched_n135_cat, matched_z_cat


def calculate_cumulative(n_cat: DataFrame, n135_cat: DataFrame, z_cat: DataFrame):
    """Works out the cumulative plots for the combined and split z selection."""

    counts = []
    non_z_counts = []
    z_counts = []
    for i in depth:
        cut_n = np.where(n_cat['MAG_CORR'] < N_LIM + i)[0]
        cut_135 = np.where(n135_cat['MAG_CORR'] < N135_LIM + i)[0]
        detected = np.intersect1d(cut_n, cut_135)
        counts.append(len(detected))
        cut_z =np.where(z_cat['MAG_CORR'] > Z_LIM)[0]
        no_z_detection = np.intersect1d(detected, cut_z)
            
        non_z_counts.append(len(no_z_detection))
        z_counts.append(len(detected) - len(no_z_detection))
        print(len(no_z_detection), len(detected), z_counts[-1])

    return np.array(counts), np.array(non_z_counts), np.array(z_counts)


if __name__ == '__main__':
    our_prefix = '../correct_stacks/N964/'
    cdfs_prefix = '../CDFS_LAGER/'

    us = prepare_data(our_prefix+'n964.cat', our_prefix+'n964_135.cat', our_prefix+'z.cat', 'candidates.txt', zero_points)
    cdfs = prepare_data(cdfs_prefix+'n964_cdfs.cat', cdfs_prefix+'n964_135_cdfs.cat', cdfs_prefix+'z_cdfs.cat', 'candidates_cdfs.txt', zero_points_cdfs)
    imacs = prepare_data('imacs_n964.cat', 'imacs_n964_135.cat', 'imacs_z.cat', 'candidates_imacs.txt', zero_points)

    counts_imacs, non_z_counts_imacs, z_counts_imacs = calculate_cumulative(imacs[0], imacs[1], imacs[2])
    counts, non_z_counts, z_counts = calculate_cumulative(us[0], us[1], us[2])
    counts_cdfs, non_z_counts_cdfs, z_counts_cdfs = calculate_cumulative(cdfs[0], cdfs[1], cdfs[2])

    #Plotting
    # Run the identify candidates scripts and update the cdfs candidates.
    #IMACS
    start_plot('Relative Narrow Band Selection', r'Counts per Area [arcsecond$^{-2}$]')
    plt.plot(depth, np.array(counts_imacs)/IMACS_AREA, label='This work', color='k', lw=3)
    plt.plot(depth, np.array(counts_cdfs)/CDFS_AREA, label='CDFS', color='r', ls=':', lw=2)
    plt.legend()
    end_plot('plots/cp_all_ic.png')

    start_plot('Relative Narrow Band Selection', r'Counts per Area [arcsecond$^{-2}$]')
    plt.plot(depth, np.array(non_z_counts_imacs)/IMACS_AREA, label='This work', color='k', lw=3)
    plt.plot(depth, np.array(non_z_counts_cdfs)/CDFS_AREA, label='CDFS', color='r', ls=':', lw=2)
    plt.legend()
    end_plot('plots/cp_nonz_ic.png')

    start_plot('Relative Narrow Band Selection', r'Counts per Area [arcsecond$^{-2}$]')
    plt.plot(depth, np.array(z_counts_imacs)/IMACS_AREA, label='This work', color='k', lw=3)
    plt.plot(depth, np.array(z_counts_cdfs)/CDFS_AREA, label='CDFS', color='r', ls=':', lw=2)
    plt.legend()
    end_plot('plots/cp_z_ic.png')


    ##DECAM
    '''start_plot('Relative Narrow Band Selection', r'Counts per Area [arcsecond$^{-2}$]')
    plt.plot(depth, np.array(counts)/DECAM_AREA, label='This work', color='k', lw=3)
    plt.plot(depth, np.array(counts_cdfs)/CDFS_AREA, label='CDFS', color='r', ls=':', lw=2)
    plt.legend()
    end_plot('plots/cp_all_dc.png')

    start_plot('Relative Narrow Band Selection', r'Counts per Area [arcsecond$^{-2}$]')
    plt.plot(depth, np.array(non_z_counts)/DECAM_AREA, label='This work', color='k', lw=3)
    plt.plot(depth, np.array(non_z_counts_cdfs)/CDFS_AREA, label='CDFS', color='r', ls=':', lw=2)
    plt.legend()
    end_plot('plots/cp_nonz_dc.png')

    start_plot('Relative Narrow Band Selection', r'Counts per Area [arcsecond$^{-2}$]')
    plt.plot(depth, np.array(z_counts)/DECAM_AREA, label='This work', color='k', lw=3)
    plt.plot(depth, np.array(z_counts_cdfs)/CDFS_AREA, label='CDFS', color='r', ls=':', lw=2)
    plt.legend()
    end_plot('plots/cp_z_dc.png')'''
