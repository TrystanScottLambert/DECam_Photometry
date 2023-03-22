""" Module to plot postage stamp across all three filters based on positon. """

from astropy.io import fits
from astropy.wcs import WCS
import pylab as plt

IMAGES = (
    '../correct_stacks/N964/i.fits',
    '../correct_stacks/N964/z.fits',
    '../correct_stacks/N964/n964.fits',
)

FITS_OBJECTS = [fits.open(image) for image in IMAGES]

PAD = 30 # pixels

def cut_postage_stamp(r_a: float, dec: float, image):
    """Cutting out a postage stamp centered on r_a, and dec (in decimal degrees)"""
    wcs = WCS(image[0].header)
    x_pix, y_pix  = wcs.world_to_pixel_values(r_a, dec)
    data = image[0].data[
        int(y_pix)-PAD:int(y_pix)+PAD, int(x_pix)-PAD:int(x_pix)+PAD]
    return data

def show_stamps(r_a: float, dec: float):
    """Cuts out postage stamps for all three images."""
    data_i = cut_postage_stamp(r_a, dec, FITS_OBJECTS[0])
    data_z = cut_postage_stamp(r_a, dec, FITS_OBJECTS[1])
    data_n964 = cut_postage_stamp(r_a, dec, FITS_OBJECTS[2])

    fig = plt.figure()
    #print(f'{r_a} {dec}')
    ax_i = fig.add_subplot(131)
    ax_i.imshow(data_i)
    plt.title('i-band')
    ax_z = fig.add_subplot(132)
    ax_z.imshow(data_z)
    plt.title('z-band')
    ax_n964 = fig.add_subplot(133)
    ax_n964.imshow(data_n964)
    plt.title('NB964-band')
    plt.show()

if __name__ == '__main__':
    RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
    DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1

    show_stamps(RA_QSO, DEC_QSO)
    