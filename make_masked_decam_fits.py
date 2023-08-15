"""
Making a masked fits file for our decam field
"""

from astropy.io import fits
from regionfy_catalog import load_region

FITS_FILE = '../correct_stacks/N964/n964.fits'
REGION_FILE = 'decam.reg'

image = fits.open(FITS_FILE)
region_decam_fov = load_region(REGION_FILE)
print('making mask')
region_decam_mask = region_decam_fov.to_mask()
print('converting to image')
region_image = region_decam_mask.to_image(image[0].data.shape)
print('mask done, writing to file')
image[0].data = region_image
image.writeto('DECAM_MASK.fits', overwrite=True)
