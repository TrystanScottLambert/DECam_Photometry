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

PAD = 30

def cut_postage_stamp(r_a: float, dec: float, image, pad = PAD):
    """Cutting out a postage stamp centered on r_a, and dec (in decimal degrees)"""
    if image[0].data is not None:
        wcs = WCS(image[0].header)
        x_pix, y_pix  = wcs.world_to_pixel_values(r_a, dec)
        data = image[0].data[
            int(y_pix)-pad:int(y_pix)+pad, int(x_pix)-pad:int(x_pix)+pad]
    else:
        wcs = WCS(image[1].header, naxis=2)
        x_pix, y_pix  = wcs.world_to_pixel_values(r_a, dec)
        data = image[1].data[
            int(y_pix)-pad:int(y_pix)+pad, int(x_pix)-pad:int(x_pix)+pad]
    return data

def cut_out_array(array, x_pix, y_pix, pad):
    """Cuts out a postage stamp of an array."""
    return array[int(y_pix)-pad:int(y_pix)+pad, int(x_pix)-pad:int(x_pix)+pad]

def cut_out_stamps(r_a: float, dec: float, fits_objects: list, **kwargs):
    """Cuts out all of the images """
    data_i = cut_postage_stamp(r_a, dec, fits_objects[0], **kwargs)
    data_z = cut_postage_stamp(r_a, dec, fits_objects[1], **kwargs)
    data_n964 = cut_postage_stamp(r_a, dec, fits_objects[2], **kwargs)
    return data_i, data_z, data_n964

def get_data_and_wcs(fits_object):
    """Returns the data and the wcs."""
    if fits_object[0].data is not None:
        hdr, data = fits_object[0].header, fits_object[0].data
    else:
        hdr, data = fits_object[1].header, fits_object[1].data
    return data, WCS(hdr)

def cut_out_mulitple_stamps(ra_array, dec_array, fits_objects, pad = 30):
    """Does the cutting out for a whole range of images."""

    all_data_i, wcs_i = get_data_and_wcs(fits_objects[0])
    all_data_z, wcs_z = get_data_and_wcs(fits_objects[1])
    all_data_n964, wcs_n964 = get_data_and_wcs(fits_objects[2])

    i_band, z_band, n964_band = [], [], []
    for i, _ in enumerate(ra_array):
        x_i, y_i = wcs_i.world_to_pixel_values(ra_array[i], dec_array[i])
        cut_out_i = cut_out_array(all_data_i, x_i, y_i, pad)
        i_band.append(cut_out_i)

        x_z, y_z = wcs_z.world_to_pixel_values(ra_array[i], dec_array[i])
        cut_out_z = cut_out_array(all_data_z, x_z, y_z, pad)
        z_band.append(cut_out_z)

        x_n964, y_n964 = wcs_n964.world_to_pixel_values(ra_array[i], dec_array[i])
        cut_out_n964 = cut_out_array(all_data_n964, x_n964, y_n964, pad)
        n964_band.append(cut_out_n964)
    return i_band, z_band, n964_band


def show_stamps(r_a: float, dec: float):
    """Cuts out postage stamps for all three images."""
    data_i = cut_postage_stamp(r_a, dec, FITS_OBJECTS[0])
    data_z = cut_postage_stamp(r_a, dec, FITS_OBJECTS[1])
    data_n964 = cut_postage_stamp(r_a, dec, FITS_OBJECTS[2])

    fig = plt.figure()
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

    i, z, n964 = cut_out_stamps(RA_QSO, DEC_QSO, FITS_OBJECTS, pad=50)
    plt.imshow(i)
    plt.show()
    plt.imshow(z)
    plt.show()
    plt.imshow(n964)
    plt.show()

    show_stamps(RA_QSO, DEC_QSO)
    