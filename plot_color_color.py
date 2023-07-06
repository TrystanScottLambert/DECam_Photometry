"""
Plottting color-color diagram using the colors set 
by eduardo and chira in 2013 and 2017.
"""
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

from sex_catalog import SExtractorCat
from zero_points import zero_points, ZeroPoint, ZeroPoints
from identify_candidates import write_region_file


def cross_match_to_sexcat(
        ra_array: np.ndarray, dec_array:np.ndarray, sex_catalog: DataFrame
        ) -> DataFrame:
    """Reduces the sexcatalog to the cross matched values"""
    sex_ra, sex_dec = np.array(
        sex_catalog['ALPHAPEAK_J2000']), np.array(sex_catalog['DELTAPEAK_J2000'])
    catalog = SkyCoord(ra = sex_ra*u.deg, dec = sex_dec*u.deg)
    candidate_coords = SkyCoord(ra = ra_array * u.deg, dec = dec_array*u.deg)
    idx, _, _ = candidate_coords.match_to_catalog_sky(catalog)
    return sex_catalog.iloc[idx]

def add_mag(catalog: DataFrame, zero_point_function: ZeroPoint) -> None:
    """Adds AB mag data to the catlog based on the zero point."""
    catalog['MAG_CORR'] = catalog['MAG_APER'] + zero_point_function.mag_correct(1)


def set_non_detections(catalog: DataFrame, three_sigma_value: float) -> None:
    """
    Sets the sextractor value of 99 to the 3sigma value, giving
    us upper and lower limits in the plot.
    """
    catalog.loc[catalog['MAG_APER'] == 99, 'MAG_CORR'] = three_sigma_value


def read_in_catalogs(
        catalog_names: list[str], zero_points_func: ZeroPoints, three_sigma_depth: list[float]
        ) -> tuple[DataFrame, DataFrame, DataFrame]:
    """
    Reads in the i, z, and narrow band catalogs returning them as dataframes.
    The lists needs to be passed in the order i, z, and narrowband.
    """
    i_cat = SExtractorCat(catalog_names[0]).catalog
    z_cat = SExtractorCat(catalog_names[1]).catalog
    n_cat = SExtractorCat(catalog_names[2]).catalog

    add_mag(i_cat, zero_points_func.i_band)
    add_mag(z_cat, zero_points_func.z_band)
    add_mag(n_cat, zero_points_func.n964_band)

    set_non_detections(i_cat, three_sigma_depth[0])
    set_non_detections(z_cat, three_sigma_depth[1])
    set_non_detections(n_cat, three_sigma_depth[2])
    return i_cat, z_cat, n_cat


if __name__ == '__main__':
    OUR_DEPTHS = [25.66, 25.58, 24.66]
    OUR_CATALOGS = [
        '../correct_stacks/N964/i.cat',
        '../correct_stacks/N964/z.cat',
        '../correct_stacks/N964/n964.cat']
    CANDIDATES_FILE = 'candidates_e.txt'
    ra, dec = np.loadtxt(CANDIDATES_FILE, unpack=True)

    our_i_cat, our_z_cat, our_n964_cat = read_in_catalogs(OUR_CATALOGS, zero_points, OUR_DEPTHS)
    cats = [our_i_cat, our_z_cat, our_n964_cat]
    candidate_cats = [cross_match_to_sexcat(ra, dec, cat) for cat in cats]
    candidate_z_nb = candidate_cats[1]['MAG_CORR'] - candidate_cats[2]['MAG_CORR']
    candidate_i_z = candidate_cats[0]['MAG_CORR'] - candidate_cats[1]['MAG_CORR']


    z_nb = our_z_cat['MAG_CORR'] - our_n964_cat['MAG_CORR']
    i_z  = our_i_cat['MAG_CORR'] - our_z_cat['MAG_CORR']
    limited_z = np.where(our_z_cat['MAG_APER'] == 99)[0]
    limited_i = np.where(our_i_cat['MAG_APER'] == 99)[0]
    limited_n964 = np.where(our_n964_cat['MAG_APER'] == 99)[0]
    only_z = np.setdiff1d(limited_z, limited_i)
    double_uncertainty = np.intersect1d(limited_i,limited_z)
    only_i = np.setdiff1d(limited_i, limited_z)
    plt.axhline(0.75, color='r', ls='--')
    plt.axvline(1, color='r', ls='--')
    plt.scatter(i_z, z_nb, color='k', s=1)
    arrow_err = (2.5/np.log(10)) / 5
    #plt.errorbar(z_nb[only_z], i_z[only_z], color='b', yerr=arrow_err, fmt='o', ms=2,elinewidth=1, lolims=True)
    #plt.errorbar(z_nb[only_i], i_z[only_i], color='g', xerr=arrow_err, fmt='o', ms=2, elinewidth=1, xlolims=True)
    #plt.errorbar(z_nb[double_uncertainty], i_z[double_uncertainty], color='r', xerr=arrow_err, yerr=arrow_err, fmt='o', ms=2, elinewidth=1, xlolims=True, lolims=True)
    plt.scatter(candidate_i_z, candidate_z_nb, marker = 'd', s=80, color='r')
    plt.show()

    write_region_file(list(our_n964_cat['ALPHAPEAK_J2000'][only_z]), list(our_n964_cat['DELTAPEAK_J2000'][only_z]), 'delete_z.reg')
