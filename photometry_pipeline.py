"""Pipeline to put DECAM data in correct format
This includes: background Subtraction and converting,
counts to flux (electrons per second)."""

from astropy.io import fits
from astropy.stats import SigmaClip
from photutils.background import Background2D, SExtractorBackground
from tqdm import tqdm

def add_pipe_method_to_name(previous_name, pipe_method):
    """Takes the name and appends a method so that we know what process has taken place."""
    prefix = previous_name.split('.fits.fz')[0]
    new_name = f'{prefix}_{pipe_method}.fits.fz'
    return new_name

def calculate_background(fits_image_data):
    """Identifies the background in the image. This would be a single tile."""
    sigma_clip = SigmaClip(sigma=3.)
    sky_estimator = SExtractorBackground()
    sky = Background2D(fits_image_data, fits_image_data.shape, filter_size = (3, 3),
                   sigma_clip = sigma_clip, bkg_estimator = sky_estimator)
    return sky

def remove_background(fits_image_data):
    """Works out the sky value and removes it from the data."""
    sky = calculate_background(fits_image_data)
    print('removing background...')
    subtracted_data = fits_image_data - sky.background
    return subtracted_data

def remove_background_from_mosaic(fits_mosaic_name):
    """Go through entire mosaic and remove the background."""
    print(f'removing sky from {fits_mosaic_name}')
    hdu = fits.open(fits_mosaic_name)
    for i in tqdm(range(1, len(hdu))): #ignore the first frame which is the header.
        hdu[i].data = remove_background(hdu[i].data)
    hdu.writeto(add_pipe_method_to_name(fits_mosaic_name, 'skysubtracted'), overwrite = True)

'''def divide_by_exposure_time(image_mosaic_name, exposure_mosaic_name):
    """Goes through a mosaic and divides by the exposure time, which results in flux."""
    print(f'Dividing {image_mosaic_name} through by exposure {exposure_mosaic_name}')
    hdu_image = fits.open(image_mosaic_name)
    hdu_exptime = fits.open(exposure_mosaic_name)
    for i in tqdm(range(1, len(hdu_image))):
        hdu_exptime[i].data[hdu_exptime[i].data == 0] = 1
        hdu_image[i].data = hdu_image[i].data / hdu_exptime[i].data
    hdu_image.writeto(add_pipe_method_to_name(image_mosaic_name, 'flux'), overwrite = True)'''

def divide_by_exposure_time(fits_mosaic_name):
    """Dividing by the exposure time"""
    hdu = fits.open(fits_mosaic_name)
    time = hdu[0].header['EXPTIME']
    for i in range(1, len(hdu)):
        hdu[i].data = hdu[i].data / time
    hdu.writeto(add_pipe_method_to_name(fits_mosaic_name, 'flux'), overwrite = True)


if __name__ == '__main__':
    FOLDER = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/'
    IMAGE_NAME = 'c4d_211021_003940_osj_i_vik1.fits.fz'
    EXP_NAME = 'c4d_211021_003940_ose_i_vik1.fits.fz'
    divide_by_exposure_time(FOLDER + 'c4d_211021_003940_osj_i_vik1_skysubtracted.fits.fz')
