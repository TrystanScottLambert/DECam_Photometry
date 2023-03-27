"""Trying to remove the cosmic rays from the images."""

from astropy.io import fits
import lacosmic

if __name__ == '__main__':
    INFILE = '../correct_stacks/N964/n964.fits'
    hdu = fits.open(INFILE)
    msk, array = lacosmic.lacosmic(hdu[0].data)
    hdu[0].data = array
    hdu.writeto('n964_remove_cosmics.fits')
