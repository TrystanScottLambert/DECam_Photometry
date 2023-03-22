"""Masking the very negative pixels of the DECam image before running through pipeline."""

import numpy as np
from astropy.io import fits


def convert_neg_to_nan(array_2d):
    """Returns a 2d array with the very negative values masked with nans."""
    neg_idx = np.where(array_2d < 0)
    array_2d[neg_idx] = np.nan
    return array_2d

def mask_neg_values(hdu_list):
    """Loops through the entire hdu list and masks the very negative values."""
    for i in range(1, len(hdu_list)): # 0th value is the headers
        hdu_list[i].data = convert_neg_to_nan(hdu_list[i].data)
    return hdu_list

def mask_fits_file(fits_name: str):
    """masks the very negative values of a given fits file."""
    hdu = fits.open(fits_name)
    hdu = mask_neg_values(hdu)
    hdu.writeto(fits_name, overwrite=True)

if __name__ == '__main__':
    INFILE_N964 = '../correct_stacks/N964/c4d_210831_050404_osj_N964_vik1.fits'
    INFILE_I = '../correct_stacks/i/c4d_211021_003940_osj_i_vik1.fits'
    INFILE_Z = '../correct_stacks/z/c4d_210831_053503_osj_z_vik1.fits'
    infiles = [INFILE_N964, INFILE_I, INFILE_Z]

    for file in infiles:
        mask_fits_file(file)
