"""
Script for working out area used for source detection. Note that this will take a while.
"""

import numpy as np
from rich.progress import track
from regions import PixCoord
from mask_sex_catalog import create_weight_mask
from regionfy_catalog import load_region

PIXELS_IN_REGION = 4.016066e8
INFILE = 'DECAM.reg'
region  = load_region(INFILE)

WEIGHTS = [
    ('../correct_stacks/N964/n964_weight.fits', 0.008),
    ('../correct_stacks/N964/i_weight.fits', 0.0008),
    ('../correct_stacks/N964/z_weight.fits', 0.001)]

msk = create_weight_mask(WEIGHTS)

x_bad, y_bad = np.where(msk==0)
pixcoords = PixCoord(x_bad, y_bad)
in_region_array = np.array([pixcoord in region for pixcoord in track(pixcoords)])
number_bad_in_region = len(np.where(in_region_array == True)[0])
print(number_bad_in_region)
total_pixels = PIXELS_IN_REGION - number_bad_in_region
print('Total number of pixels are', total_pixels)
