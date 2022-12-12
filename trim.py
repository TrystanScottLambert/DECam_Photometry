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

def cutout_position(hdu, position: SkyCoord, size: float, wcs: WCS):
    """Cutsout the position"""
    cutout = Cutout2D(hdu[0].data, position, size, wcs=wcs)
    return cutout

def make_new_fits_file(hdu, cutout_object: Cutout2D, outfile: str):
    """Writes a new fits object based on the cutout data."""
    new_hdu = fits.PrimaryHDU()
    new_hdu.data = cutout_object.data
    new_hdu.header = hdu[0].header
    new_hdu.header.update(cutout.wcs.to_header())
    new_hdu.writeto(outfile, overwrite=True)


if __name__ == '__main__':
    INFILE_I = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/i_band_coadd.fits'
    INFILE_Z = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/z/z_band_coadd.fits'
    INFILE_N964 = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/N964/N964_band_coadd.fits'
    SIZE = 2*u.deg

    hdu, wcs = read_in_fits(INFILE_N964)
    position = SkyCoord(hdu[0].header['CRVAL1']*u.deg, hdu[0].header['CRVAL2']*u.deg)
    cutout = cutout_position(hdu, position, SIZE, wcs)
    make_new_fits_file(hdu, cutout, 'test_delete_N964.fits')

    hdu_i, wcs_i = read_in_fits(INFILE_I)
    cutout_i = cutout_position(hdu_i, position, SIZE, wcs_i)
    make_new_fits_file(hdu_i, cutout_i, 'test_delete_i.fits')

    hdu_z, wcs_z = read_in_fits(INFILE_Z)
    cutout_z = cutout_position(hdu_z, position, SIZE, wcs_z)
    make_new_fits_file(hdu_z, cutout_z, 'test_delete_z.fits')
