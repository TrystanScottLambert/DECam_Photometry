"""
This script trims FITS images for sextractor comparative photometry analysis
by creating a new FITS object based on a
cutout of the original data.
"""

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


def read_in_fits(infile: str):
    """
    Reads in the FITS file from the given input file path.

    Args:
        infile (str): Input file path.

    Returns:
        Tuple containing the HDUList object and the WCS object extracted from the header.
    """
    hdu = fits.open(infile)
    wcs = WCS(hdu[0].header)
    return hdu, wcs

def make_new_fits_file(hdu:fits.hdu.hdulist.HDUList, cutout_object: Cutout2D, outfile: str):
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
    new_hdu.header = hdu[0].header
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
    cutout = Cutout2D(hdu[0].data, position, size, wcs=wcs)
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
    center = SkyCoord(hdu[0].header['CRVAL1'] * u.deg, hdu[0].header['CRVAL2']*u.deg)
    return center

if __name__ == '__main__':
    INFILE_I = '../correct_stacks/i/i_band_coadd.fits'
    INFILE_I_WEIGHT = '../correct_stacks/i/i_band_weight_coadd.fits'

    INFILE_Z = '../correct_stacks/z/z_band_coadd.fits'
    INFILE_Z_WEIGHT = '../correct_stacks/z/z_band_weight_coadd.fits'

    INFILE_N964 = '../correct_stacks/N964/N964_band_coadd.fits'
    INFILE_N964_WEIGHT = '../correct_stacks/N964/N964_band_weight_coadd.fits'
    SIZE = 2 * u.deg

    n964_position = get_center_of_image(INFILE_N964)

    cut_fits_image(INFILE_I, n964_position, SIZE, '../correct_stacks/N964/i.fits')
    cut_fits_image(INFILE_I_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/i_weight.fits')

    cut_fits_image(INFILE_Z, n964_position, SIZE, '../correct_stacks/N964/z.fits')
    cut_fits_image(INFILE_Z_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/z_weight.fits')

    cut_fits_image(INFILE_N964, n964_position, SIZE, '../correct_stacks/N964/n964.fits')
    cut_fits_image(INFILE_N964_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/n964_weight.fits')
