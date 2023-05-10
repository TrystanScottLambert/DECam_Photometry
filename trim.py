"""
This script trims FITS images for sextractor comparative photometry analysis
by creating a new FITS object based on a
cutout of the original data.
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


def read_in_fits(infile: str) -> tuple[fits.HDUList, WCS]:
    """
    Reads in the FITS file from the given input file path.

    Args:
        infile (str): Input file path.

    Returns:
        Tuple containing the HDUList object and the WCS object extracted from the header.
    """
    hdu = fits.open(infile)
    wcs = WCS(hdu[1].header)
    return hdu, wcs

def make_new_fits_file(hdu:fits.hdu.hdulist.HDUList, cutout_object: Cutout2D, outfile: str) -> None:
    """
    Writes a new FITS file based on the cutout data and header information from the original file.

    Args:
        hdu (HDUList): Original FITS file's HDUList object.
        cutout_object (Cutout2D): Cutout object created using
        original FITS data and WCS information.
        outfile (str): Output file path.
    """
    new_hdu = fits.PrimaryHDU()
    new_hdu.data = cutout_object.data
    new_hdu.header = hdu[1].header
    new_hdu.header.update(cutout_object.wcs.to_header())
    new_hdu.writeto(outfile, overwrite=True)

def cut_fits_image(fits_image: str, position: SkyCoord, size: float, outfile: str) -> None:
    """
    Cuts a FITS image based on a central point and a box size, and writes the new FITS file.

    Args:
        fits_image (str): Input FITS file path.
        position (SkyCoord): Central position of the box.
        size (float): Size of the box in arcseconds.
        outfile (str): Output file path.
    """
    hdu, wcs = read_in_fits(fits_image)
    cutout = Cutout2D(hdu[1].data, position, size, wcs=wcs)
    make_new_fits_file(hdu, cutout, outfile)

def get_center_of_image(fits_name: str) -> SkyCoord:
    """
    Given a FITS file name, opens and extracts the center coordinate of the image.
    
    Parameters:
    -----------
    fits_name : str
        The name of the FITS file to open and extract the center coordinate.
        
    Returns:
    --------
    center : `astropy.coordinates.SkyCoord`
        The center coordinate of the image in `astropy.coordinates.SkyCoord` format.
    """
    hdu = fits.open(fits_name)
    data = hdu[1].data
    y_pix, x_pix = data.shape[0]/2, data.shape[1]/2
    wcs = WCS(hdu[1].header)
    ra, dec = wcs.pixel_to_world_values(x_pix, y_pix)
    center = SkyCoord(ra * u.deg, dec * u.deg)
    return center

if __name__ == '__main__':
    INFILE_I = '../correct_stacks/i/c4d_211021_003551_osj_i_vik2.fits.fz'
    INFILE_I_WEIGHT = '../correct_stacks/i/c4d_211021_003551_osw_i_vik2.fits.fz'

    INFILE_Z = '../correct_stacks/z/c4d_210831_053503_osj_z_vik2.fits.fz'
    INFILE_Z_WEIGHT = '../correct_stacks/z/c4d_210831_053503_osw_z_vik2.fits.fz'

    INFILE_N964 = '../correct_stacks/N964/c4d_210831_050404_osj_N964_vik2.fits.fz'
    INFILE_N964_WEIGHT = '../correct_stacks/N964/c4d_210831_050404_osw_N964_vik2.fits.fz'
    SIZE = 1.5 * u.deg

    n964_position = get_center_of_image(INFILE_N964)

    cut_fits_image(INFILE_I, n964_position, SIZE, '../correct_stacks/N964/i.fits')
    cut_fits_image(INFILE_I_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/i_weight.fits')

    cut_fits_image(INFILE_Z, n964_position, SIZE, '../correct_stacks/N964/z.fits')
    cut_fits_image(INFILE_Z_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/z_weight.fits')

    cut_fits_image(INFILE_N964, n964_position, SIZE, '../correct_stacks/N964/n964.fits')
    cut_fits_image(INFILE_N964_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/n964_weight.fits')
