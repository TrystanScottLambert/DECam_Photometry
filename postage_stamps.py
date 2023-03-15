""" Module to plot postage stamp across all three filters based on positon. """

from astropy.io import fits
from astropy.wcs import WCS
import pylab as plt
import astropy
IMAGES = (
    '../correct_stacks/N964/i.fits',
    '../correct_stacks/N964/z.fits',
    '../correct_stacks/N964/n964.fits',
)

FITS_OBJECTS = [fits.open(image) for image in IMAGES]

PAD = 20 # pixels

def cut_postage_stamp(ra: float, dec: float, image):
    """Cutting out a postage stamp centered on ra, and dec (in decimal degrees)"""
    wcs = WCS(image[0].header)
    x_pix, y_pix  = wcs.world_to_pixel_values(ra, dec)
    data = image[0].data[
        int(y_pix)-PAD:int(y_pix)+PAD, int(x_pix)-PAD:int(x_pix)+PAD]
    return data

def show_stamps(ra: float, dec: float):
    """Cuts out postage stamps for all three images."""
    data_i = cut_postage_stamp(ra, dec, FITS_OBJECTS[0])
    data_z = cut_postage_stamp(ra, dec, FITS_OBJECTS[1])
    data_n964 = cut_postage_stamp(ra, dec, FITS_OBJECTS[2])

    plt.subplot(131)
    plt.imshow(data_i)
    plt.title('i-band')
    plt.subplot(132)
    plt.imshow(data_z)
    plt.title('z-band')
    plt.subplot(133)
    plt.imshow(data_n964)
    plt.title('NB964-band')
    plt.show()

if __name__ == '__main__':
    RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
    DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1

    show_stamps(RA_QSO, DEC_QSO)
    