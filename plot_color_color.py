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
from zero_points_cdfs import zero_points_cdfs
from identify_candidates import write_region_file
from snr_fit import calculate_snr


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

def read_in_catalogs(
        catalog_names: list[str], zero_points_func: ZeroPoints) -> tuple[DataFrame, DataFrame, DataFrame]:
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
    return i_cat, z_cat, n_cat


if __name__ == '__main__':
    I_BAND_2_SIGMA = 26.64
    Z_BAND_2_SIGMA = 26.58
    N_BAND_2_SIGMA = 25.69
    OUR_CATALOGS = [
        '../correct_stacks/N964/i.cat',
        '../correct_stacks/N964/z.cat',
        '../correct_stacks/N964/n964.cat']
    CANDIDATES_FILE = 'candidates_e.txt'

    #I_BAND_2_SIGMA = 28.10
    #Z_BAND_2_SIGMA = 27.73
    #N_BAND_2_SIGMA = 25.82
    #OUR_CATALOGS = [
    #    '../CDFS_LAGER/i_cdfs.cat',
    #    '../CDFS_LAGER/z_cdfs.cat',
    #    '../CDFS_LAGER/n964_cdfs.cat']
    #CANDIDATES_FILE = 'candidates_e.txt'



    ra, dec = np.loadtxt(CANDIDATES_FILE, unpack=True)

    our_i_cat, our_z_cat, our_n964_cat = read_in_catalogs(OUR_CATALOGS, zero_points)
    z_mag, z_err = our_z_cat['MAG_CORR'], our_z_cat['MAGERR_APER']
    i_mag, i_err = our_i_cat['MAG_CORR'], our_i_cat['MAGERR_APER']
    z_snr, i_snr = calculate_snr(z_err), calculate_snr(i_err)

    z_cut = np.where(our_z_cat['MAG_APER'] == 99)[0]
    i_cut = np.where(our_i_cat['MAG_APER'] == 99)[0]

    artifacts = np.intersect1d(z_cut, i_cut) # these tend to be bright stars at the saturation limit on narrow band

    n_mag, n_err = our_n964_cat['MAG_CORR'], our_n964_cat['MAGERR_APER']
    n_snr = calculate_snr(n_err)
    ra_all, dec_all = np.array(our_n964_cat['ALPHAPEAK_J2000']), np.array(our_n964_cat['DELTAPEAK_J2000'])
    
    plt.scatter(ra_all, dec_all, c = np.array(z_mag))
    plt.show()

    ra, dec = np.array(our_n964_cat['ALPHAPEAK_J2000']), np.array(our_n964_cat['DELTAPEAK_J2000'])
    cats = [our_i_cat, our_z_cat, our_n964_cat]

    color_zn = z_mag - n_mag
    color_iz = i_mag - z_mag

    cut1 = np.where(color_iz > 4)[0]
    cut2 = np.where(color_zn < 0.35)[0]
    more_artifacts = np.intersect1d(cut1, cut2) # tend to be stars saturating in i-band
    artifacts = np.union1d(artifacts, more_artifacts)
    good = np.setdiff1d(np.arange(len(i_mag)), artifacts)

    n_mag, n_err, z_mag, z_err, i_mag, i_err = n_mag[good], n_err[good], z_mag[good], z_err[good], i_mag[good], i_err[good]
    n_snr, z_snr, i_snr = n_snr[good], z_snr[good], i_snr[good]
    ra, dec = ra[good], dec[good]


    z_mag[z_snr < 2] = Z_BAND_2_SIGMA
    i_mag[i_snr < 2] = I_BAND_2_SIGMA
    print('non detections: ', len(np.where(z_snr<2)[0])/len(z_snr))
    color_zn = z_mag - n_mag
    color_iz = i_mag - z_mag

    cut1 = np.where(color_zn > 0.78)[0]
    cut2 = np.where(color_iz > 1)[0]

    pass_color_cut = np.intersect1d(cut1, cut2)
    plt.scatter(color_iz, color_zn, s=1)
    #plt.scatter(np.array(color_iz)[pass_color_cut], np.array(color_zn)[pass_color_cut])
    plt.axhline(0.78, color='r', ls='--')
    plt.axvline(1, color='r', ls='--')
    plt.xlim(-1, 2.5)
    plt.ylim(-2, 4)
    plt.show()

    z_mag_passes_cut = np.array(z_mag)[pass_color_cut]
    z_err_passes_cut = np.array(z_err)[pass_color_cut]
    z_snr_passes_cut = np.array(z_snr)[pass_color_cut]

    i_mag_passes_cut = np.array(i_mag)[pass_color_cut]
    i_err_passes_cut = np.array(i_err)[pass_color_cut]
    i_snr_passes_cut = np.array(i_snr)[pass_color_cut]

    n_mag_passes_cut = np.array(n_mag)[pass_color_cut]
    n_err_passes_cut = np.array(n_err)[pass_color_cut]
    n_snr_passes_cut = np.array(n_snr)[pass_color_cut]

    no_i_detection = np.where(i_snr_passes_cut < 2.)[0]

    z_no_i_mag, z_no_i_err = z_mag_passes_cut[no_i_detection], z_err_passes_cut[no_i_detection]
    n_no_i_mag, n_no_i_err = n_mag_passes_cut[no_i_detection], n_err_passes_cut[no_i_detection]

    color = np.abs(z_no_i_mag - n_no_i_mag)
    sig = 2.5 * np.hypot(n_no_i_err, z_no_i_err)

    final_cut = np.where(color > sig)[0]
    ra, dec = ra[pass_color_cut][no_i_detection][final_cut], dec[pass_color_cut][no_i_detection][final_cut]
    write_region_file(ra, dec, 'vis_delete.reg',size=10)

    #candidate_cats = [cross_match_to_sexcat(ra, dec, cat) for cat in cats]
    #can_z = candidate_cats[1]['MAG_CORR']
    #can_z_err = candidate_cats[1]['MAGERR_APER']
    #can_i = candidate_cats[0]['MAG_CORR']
    #can_i_err = candidate_cats[0]['MAGERR_APER']
    #can_n = candidate_cats[2]['MAG_CORR']
    #can_n_err = candidate_cats[2]['MAGERR_APER']

    #can_i[can_i>26.68] = 26.68

    #candidate_z_nb = can_z - can_n
    #candidate_i_z = can_i - can_z

    #candidate_z_nb_err = np.hypot(can_z_err, can_n_err)
    #candidate_i_z_err = np.hypot(can_z_err, can_i_err)

    #plt.axhline(0.78, color='r', ls='--')
    #plt.axvline(1, color='r', ls='--')
    #plt.scatter(color_iz[good], color_zn[good], s=1, alpha=0.3, color='k')
    #plt.errorbar(candidate_i_z, candidate_z_nb, color='r', yerr=candidate_z_nb_err, fmt='o')
    #plt.xlim(-0.7, 4.1)
    #plt.ylim(-0.7, 3.6)
    #plt.show()
