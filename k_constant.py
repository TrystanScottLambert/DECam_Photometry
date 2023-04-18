"""
Module for the K-constant which is used for determining the magnitude correction (zpt)
based on the size of the aperture which sextractor uses.
"""

import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from scipy.integrate import quad

def gaussian(x_array: np.ndarray, mean: float, fwhm: float):
    """
    Defintion of a gaussian centered on mu with a given fwhm.
    We choose to define the gaussian with a fwhm max because this 
    is equivalent to the seeing of each image.
    """
    sig = fwhm * gaussian_fwhm_to_sigma
    print(sig)
    return np.exp(-np.power(x_array - mean, 2.) / (2 * np.power(sig, 2.)))

infinite_integral = quad(gaussian, 0, np.inf, args=(0, 1.32))
definite_integral = quad(gaussian, 0, 1.35, args=(0, 1.32))

k = infinite_integral[0]/definite_integral[0]