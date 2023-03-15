"""Module to cnonvert .cat files into .reg files."""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

def read_in_positions(file_name: str):
    """Reads in the first and second column which are ra and dec respectively"""
    ra, dec = np.loadtxt(file_name, usecols=(0,1), unpack=True)
    return ra, dec

def replace_letter_with_colons(string: str):
    """30d45m64s --> 30:45:64"""
    letters = ('h','d','m')
    for letter in letters:
        string = string.replace(letter, ':')
    return string.replace('s','')

def convert_decimal_degrees_into_celestial(ra: float, dec: float):
    """Takes floats and returns them in the hh:mm:ss and dd:mm:ss format."""
    c = SkyCoord(ra = ra*u.deg, dec = dec*u.deg)
    return replace_letter_with_colons(c.to_string('hmsdms'))


def convert_cat_to_reg(sextractor_catalog: str, outfile: str = 'sex_cat.reg'):
    """takes a .cat file and returns a .reg file"""
    ra, dec = read_in_positions(sextractor_catalog)
    positions = [convert_decimal_degrees_into_celestial(ra[i], dec[i]) for i in range(len(ra))]
    f = open(outfile, 'w')
    for pos in positions:
        f.write(f'circle {pos} 2.5" \n')
    f.close()


if __name__ == '__main__':
    INFILE = '../correct_stacks/N964/n964.cat'
    convert_cat_to_reg(INFILE)
