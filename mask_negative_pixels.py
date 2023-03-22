"""
This module provides functions to mask very negative pixel values in DECam FITS images.

The module contains two functions:

    - mask_very_negative_pixels: Replace very negative pixel values with NaNs in a 2D NumPy array.
    - mask_fits_file: Mask the very negative pixel values in a DECam FITS file by calling
      mask_very_negative_pixels on each HDU data array.

Usage:
------
To mask the very negative pixel values in a set of DECam FITS files,
provide a list of filepaths to mask_fits_file.

Example:
--------
Suppose you have a set of DECam FITS files:

    - image1.fits
    - image2.fits
    - image3.fits

To mask the very negative pixel values in each file, call mask_fits_file with a list of filepaths:

    from mask_neg_pixels import mask_fits_file

    input_files = ['image1.fits', 'image2.fits', 'image3.fits']
    for file_path in input_files:
        mask_fits_file(file_path)

This will create new files with the same names as the input files,
with the very negative pixel values masked with NaNs.
"""

import numpy as np
from astropy.io import fits

def mask_very_negative_pixels(data: np.ndarray) -> np.ndarray:
    """
    Replace very negative pixel values with NaNs.

    Parameters
    ----------
    data : numpy.ndarray
        2D array of pixel values.

    Returns
    -------
    numpy.ndarray
        2D array with very negative values masked with NaNs.
    """
    masked_data = np.where(data < 0, np.nan, data)
    return masked_data

def mask_fits_file(filepath: str) -> None:
    """
    Mask the very negative pixel values in a FITS file.

    Parameters
    ----------
    filepath : str
        Filepath of the FITS file to mask.

    Returns
    -------
    None
    """
    with fits.open(filepath) as hdulist:
        for hdu in hdulist[1:]:
            hdu.data = mask_very_negative_pixels(hdu.data)
        hdulist.writeto(filepath, overwrite=True)

if __name__ == '__main__':
    input_files = [
        '../correct_stacks/N964/c4d_210831_050404_osj_N964_vik1.fits',
        '../correct_stacks/i/c4d_211021_003940_osj_i_vik1.fits',
        '../correct_stacks/z/c4d_210831_053503_osj_z_vik1.fits'
    ]

    for file_path in input_files:
        mask_fits_file(file_path)
