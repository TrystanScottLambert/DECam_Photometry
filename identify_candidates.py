"""Applying the selection criteria and identifying LAE galaxy candidates.
The criteria are the same used in Hu et. al., (2019):

1. S/N_2" > 5 and S/N_1.35" > 5 for the N964 filter. 
2. S/N_2" < 3 for the i band filter.
3a. z-N964 > 1.9 and S/N_2" > 3 for the z filter or
3b. S/N_2" < 3 for the z filter. """

import numpy as np
import convert_sexcat_into_region as con


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

def write_region_file(ra_array: np.array, dec_array: np.array) -> None:
    """Writes a region file with decimal ra and dec arrays."""
    positions = [con.convert_decimal_degrees_into_celestial(ra_array[i], dec_array[i]) \
                  for i in range(len(ra_array))]
    file = open('candidates.reg', 'w', encoding='utf8')
    for pos in positions:
        file.write(f'circle {pos} 2.0" \n')
    file.close()

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

    # Write the candidates as a region file
    ra, dec = np.loadtxt(INFILE_N964_2, usecols=(0,1), unpack=True)
    ra, dec = ra[final_cut], dec[final_cut]
    write_region_file(ra, dec)
