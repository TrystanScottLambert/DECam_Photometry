"""Module to cnonvert .cat files into .reg files. Note you can also you the aperture_check.fits"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

def read_in_positions(file_name: str):
    """Reads in the first and second column which are r_a and dec respectively"""
    r_a, dec = np.loadtxt(file_name, usecols=(0,1), unpack=True)
    return r_a, dec

def replace_letter_with_colons(string: str):
    """30d45m64s --> 30:45:64"""
    letters = ('h','d','m')
    for letter in letters:
        string = string.replace(letter, ':')
    return string.replace('s','')

def convert_decimal_degrees_into_celestial(r_a: float, dec: float):
    """Takes floats and returns them in the hh:mm:ss and dd:mm:ss format."""
    c = SkyCoord(r_a = r_a*u.deg, dec = dec*u.deg)
    return replace_letter_with_colons(c.to_string('hmsdms'))


def convert_cat_to_reg(sextractor_catalog: str, outfile: str = 'sex_cat.reg'):
    """takes a .cat file and returns a .reg file"""
    r_a, dec = read_in_positions(sextractor_catalog)
    positions = [convert_decimal_degrees_into_celestial(r_a[i], dec[i]) for i in range(len(r_a))]
    file = open(outfile, 'w', encoding='utf8')
    for pos in positions:
        file.write(f'circle {pos} 2.5" \n')
    file.close()


if __name__ == '__main__':
    INFILE = '../correct_stacks/N964/n964.cat'
    convert_cat_to_reg(INFILE)
