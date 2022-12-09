""""Sort directory of data into SCIENCE BIAS AND FLAT fields.
This is being used specifically to organize the IMACS data for THELI"""

import os
from astropy.io import fits 

OPTIONS = [
    'Bias', 
    'Flat', 
    'Science',
]

def get_fits_type(fits_file: str) -> str:
    """Reads the header and returns the fits file type."""
    hdu = fits.open(fits_file)
    object_type = hdu[0].header['OBJECT']
    if object_type in OPTIONS:
        pass