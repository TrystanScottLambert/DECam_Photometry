"""Trying to remove the cosmic rays from the images."""

from astropy.io import fits
import astroscrappy as asc

if __name__ == '__main__':
    INFILE = '../correct_stacks/N964/c4d_210831_050404_osj_N964_vik1.fits'
    hdul = fits.open(INFILE)
    for hdu in hdul[1:]:
        msk, array = asc.detect_cosmics(hdu.data)
        hdu.data = array
    hdul.writeto('n964_remove_cosmics.fits')
