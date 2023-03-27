"""Applying the selection criteria and identifying LAE galaxy candidates.
The criteria are the same used in Hu et. al., (2019):

1. S/N_2" > 5 and S/N_1.35" > 5 for the N964 filter. 
2. S/N_2" < 3 for the i band filter.
3a. z-N964 > 1.9 and S/N_2" > 3 for the z filter or
3b. S/N_2" < 3 for the z filter. """

import warnings
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

import convert_sexcat_into_region as con
import postage_stamps as ps
from gui import start_gui

RED_LIST_NAME = 'candidates_red_list.txt'

def calculate_snr(mag_err: float) -> float:
    """Converts the magnitude error into a snr value."""
    return (2.5/np.log(10))/mag_err

def find_values(value: int, array: np.array, function: str = 'Greater') -> np.array:
    """Allows user to identify values in an array 'Greater' than or 'Less' than a given value."""
    if function == 'Greater':
        idx = np.where(array > value)[0]
    elif function == 'Less':
        idx = np.where(array < value)[0]
    else:
        print('Function must be "Greater" or "Less"')
    return idx

def read_all(catalog_name: str) -> tuple[np.array, np.array, np.array]:
    """Reads in the magnitudes, errors, and works out the snr."""
    mag, err = np.loadtxt(catalog_name, usecols=(4, 5), unpack=True)
    snr = calculate_snr(err)
    return mag, err, snr

def write_region_file(ra_array:np.array, dec_array:np.array, outfile:str, size:float = 2.) -> None:
    """Writes a region file with decimal ra and dec arrays."""
    positions = [con.convert_decimal_degrees_into_celestial(ra_array[i], dec_array[i]) \
                  for i in range(len(ra_array))]
    file = open(outfile, 'w', encoding='utf8')
    for pos in positions:
        file.write(f'circle {pos} {size}" # width=4\n')
    file.close()

def update_candidate_red_list(ra_array: np.ndarray, dec_array: np.ndarray):
    """Updates the red list of candidates which are banned from processing. (obvious artificats)"""
    with open(RED_LIST_NAME,'a+', encoding='utf8') as file:
        for i, _ in enumerate(ra_array):
            file.write(f'{ra_array[i]} {dec_array[i]} \n')

def get_red_values():
    """Reads in the ra and dec of the red sources. (Not allowed to be used)"""
    try:
        ra_rejected, dec_rejected = np.loadtxt(RED_LIST_NAME, unpack=True)
    except (OSError, ValueError):
        warnings.warn(f'{RED_LIST_NAME} not found or empty. Assuming there are no rejects.')
        ra_rejected, dec_rejected = np.array([]), np.array([])
    return ra_rejected, dec_rejected

def remove_bad_values(ra_array, dec_array):
    """Removes the previously rejected ra and dec values from the current candidates."""
    ra_bad, dec_bad = get_red_values()

    catalog = SkyCoord(ra = ra_array*u.deg, dec = dec_array*u.deg)
    c_bad = SkyCoord(ra = ra_bad *u.deg, dec = dec_bad * u.deg)
    idx, d2d, _ = c_bad.match_to_catalog_sky(catalog)

    msk = d2d < 1*u.arcsec
    if len(msk) == 1:  # edge case of n=1. Then value isn't read as an array but as a float.
        idx = np.array([idx])
    idx_bad = idx[msk]
    idx_good = np.setdiff1d(np.arange(len(ra_array)), idx_bad)

    return ra_array[idx_good], dec_array[idx_good]


if __name__ == '__main__':
    INFILE_N964_135 = '../correct_stacks/N964/n964_135.cat'
    INFILE_N964_2 = '../correct_stacks/N964/n964.cat'
    INFILE_I_2 = '../correct_stacks/N964/i.cat'
    INFILE_Z_2 = '../correct_stacks/N964/z.cat'

    #1. S/N_2" > 5 and S/N_1.35" > 5 for the N964 filter.
    mag_n964_2, mag_err_n964_2, snr_n964_2 = read_all(INFILE_N964_2)
    mag_n964_1, mag_err_n964_1, snr_n964_1 = read_all(INFILE_N964_135)

    first_cut_1 = find_values(5, snr_n964_1)
    first_cut_2 = find_values(5, snr_n964_2)
    first_cut = np.intersect1d(first_cut_1, first_cut_2)

    #2. S/N_2" < 3 for the i band filter.
    mag_i_2, mag_err_i_2, snr_i_2 = read_all(INFILE_I_2)
    second_cut = find_values(3, snr_i_2, 'Less')

    #3a.  z-N964 > 1.9 and S/N_2" > 3 for the z filter
    mag_z_2, mag_err_z_2, snr_z_2 = read_all(INFILE_Z_2)
    color = mag_z_2 - mag_n964_2

    third_cut_a_1 = find_values(1.9, color)
    third_cut_a_2 = find_values(3, snr_z_2)
    third_cut_a = np.intersect1d(third_cut_a_1, third_cut_a_2)

    #3.b  S/N_2" < 3 for the z filter.
    third_cut_b = find_values(3, snr_z_2, 'Less')

    # Find final candiates
    top_cuts = np.intersect1d(first_cut, second_cut)
    final_cut_a = np.intersect1d(top_cuts, third_cut_a)
    final_cut_b = np.intersect1d(top_cuts, third_cut_b)

    final_cut = np.union1d(final_cut_a, final_cut_b)
    print(f'Final candidate count is: {len(final_cut)}')

    # Visually inspecting the remaining candidates
    ra, dec = np.loadtxt(INFILE_N964_2, usecols=(0,1), unpack=True)
    ra, dec = ra[final_cut], dec[final_cut]
    ra, dec = remove_bad_values(ra, dec)
    print(f'After removing previous rejects, count is: {len(ra)}')


    i_bands, z_bands, n_bands = [], [], []
    for i, _ in enumerate(ra):
        i_filter, z_filter, n964_filter = ps.cut_out_stamps(ra[i], dec[i], pad=60)
        i_bands.append(i_filter)
        z_bands.append(z_filter)
        n_bands.append(n964_filter)

    artifacts, candidates = start_gui(i_bands, z_bands, n_bands)
    ra_rejects = ra[artifacts]
    dec_rejects = dec[artifacts]
    update_candidate_red_list(ra_rejects, dec_rejects)

    write_region_file(ra[candidates], dec[candidates], 'candidates.reg', size=8)
    with open('candidates.txt', 'w', encoding='utf8') as file:
        file.write('# RA DEC \n')
        for i, _ in enumerate(candidates):
            file.write(f'{ra[candidates][i]} {dec[candidates][i]} \n')
