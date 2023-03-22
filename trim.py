"""Trim fits images specifically for sextractor do to comparative photometry on same size image."""

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


def read_in_fits(infile: str):
    """Reads in the fits file."""
    hdu = fits.open(infile)
    wcs = WCS(hdu[0].header)
    return hdu, wcs

def make_new_fits_file(hdu, cutout_object: Cutout2D, outfile: str):
    """Writes a new fits object based on the cutout data."""
    new_hdu = fits.PrimaryHDU()
    new_hdu.data = cutout_object.data
    new_hdu.header = hdu[0].header
    new_hdu.header.update(cutout_object.wcs.to_header())
    new_hdu.writeto(outfile, overwrite=True)

def cut_fits_image(fits_image: str, position: SkyCoord, size: float, outfile: str) -> None:
    """Cuts a fits image based on a central point in a box shape defined by size."""
    hdu, wcs = read_in_fits(fits_image)
    cutout = Cutout2D(hdu[0].data, position, size, wcs=wcs)
    make_new_fits_file(hdu, cutout, outfile)


if __name__ == '__main__':
    INFILE_I = '../correct_stacks/i/i_band_coadd.fits'
    INFILE_I_WEIGHT = '../correct_stacks/i/i_band_weight_coadd.fits'

    INFILE_Z = '../correct_stacks/z/z_band_coadd.fits'
    INFILE_Z_WEIGHT = '../correct_stacks/z/z_band_weight_coadd.fits'

    INFILE_N964 = '../correct_stacks/N964/N964_band_coadd.fits'
    INFILE_N964_WEIGHT = '../correct_stacks/N964/N964_band_weight_coadd.fits'
    SIZE = 2 * u.deg

    hdu, wcs = read_in_fits(INFILE_N964)
    n964_position = SkyCoord(hdu[0].header['CRVAL1']*u.deg, hdu[0].header['CRVAL2']*u.deg)

    cut_fits_image(INFILE_I, n964_position, SIZE, '../correct_stacks/N964/i.fits')
    cut_fits_image(INFILE_I_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/i_weight.fits')

    cut_fits_image(INFILE_Z, n964_position, SIZE, '../correct_stacks/N964/z.fits')
    cut_fits_image(INFILE_Z_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/z_weight.fits')

    cut_fits_image(INFILE_N964, n964_position, SIZE, '../correct_stacks/N964/n964.fits')
    cut_fits_image(INFILE_N964_WEIGHT, n964_position, SIZE, '../correct_stacks/N964/n964_weight.fits')
