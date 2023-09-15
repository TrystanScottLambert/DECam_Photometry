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
    OUR_DEPTHS = [26.23, 25.58, 24.66]
    I_BAND_3_SIGMA = 26.23
    I_BAND_2_SIGMA = 26.68
    Z_BAND_3_SIGMA = 26.16
    OUR_CATALOGS = [
        '../correct_stacks/N964/i.cat',
        '../correct_stacks/N964/z.cat',
        '../correct_stacks/N964/n964.cat']
    CANDIDATES_FILE = 'candidates_e.txt'
    ra, dec = np.loadtxt(CANDIDATES_FILE, unpack=True)

    our_i_cat, our_z_cat, our_n964_cat = read_in_catalogs(OUR_CATALOGS, zero_points, OUR_DEPTHS)
    z_mag, z_err = our_z_cat['MAG_CORR'], our_z_cat['MAGERR_APER']
    z_mag[our_z_cat['MAG_APER'] == 99] = Z_BAND_3_SIGMA
    z_err[our_z_cat['MAG_APER'] == 99] = 0.1

    i_mag, i_err = our_i_cat['MAG_CORR'], our_i_cat['MAGERR_APER']
    i_mag[our_i_cat['MAG_APER'] == 99] = I_BAND_3_SIGMA
    i_err[our_i_cat['MAG_APER'] == 99] = 0.1

    z_cut = np.where(our_z_cat['MAG_APER'] == 99)[0]
    i_cut = np.where(our_i_cat['MAG_APER'] == 99)[0]

    artifacts = np.intersect1d(z_cut, i_cut) # these tend to be bright stars at the saturation limit on narrow band

    n_mag, n_err = our_n964_cat['MAG_CORR'], our_n964_cat['MAGERR_APER']
    cats = [our_i_cat, our_z_cat, our_n964_cat]

    color_zn = z_mag - n_mag
    color_iz = i_mag - z_mag

    cut1 = np.where(color_iz > 4)[0]
    cut2 = np.where(color_zn < 0.35)[0]
    more_artifacts = np.intersect1d(cut1, cut2) # tend to be stars saturating in i-band
    artifacts = np.union1d(artifacts, more_artifacts)
    good = np.setdiff1d(np.arange(len(i_mag)), artifacts)

    candidate_cats = [cross_match_to_sexcat(ra, dec, cat) for cat in cats]
    can_z = candidate_cats[1]['MAG_CORR']
    can_z_err = candidate_cats[1]['MAGERR_APER']
    can_i = candidate_cats[0]['MAG_CORR']
    can_i_err = candidate_cats[0]['MAGERR_APER']
    can_n = candidate_cats[2]['MAG_CORR']
    can_n_err = candidate_cats[2]['MAGERR_APER']

    can_i[can_i>26.68] = 26.68

    candidate_z_nb = can_z - can_n
    candidate_i_z = can_i - can_z

    candidate_z_nb_err = np.hypot(can_z_err, can_n_err)
    candidate_i_z_err = np.hypot(can_z_err, can_i_err)

    plt.axhline(0.78, color='r', ls='--')
    plt.axvline(1, color='r', ls='--')
    plt.scatter(color_iz[good], color_zn[good], s=1, alpha=0.3, color='k')
    plt.errorbar(candidate_i_z, candidate_z_nb, color='r', yerr=candidate_z_nb_err, fmt='o')
    plt.xlim(-0.7, 4.1)
    plt.ylim(-0.7, 3.6)
    plt.show()
