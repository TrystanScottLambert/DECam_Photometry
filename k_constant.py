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
    return np.exp(-np.power(x_array - mean, 2.) / (2 * np.power(sig, 2.)))

def calculate_k_constant(aperture_radius: float, seeing: float):
    """k(r). Make sure the aperture radius and seeing are in the same units."""
    infinite_integral = quad(gaussian, 0, np.inf, args=(0, seeing))
    definite_integral = quad(gaussian, 0, aperture_radius, args=(0, seeing))
    constant = infinite_integral[0]/definite_integral[0]
    return constant

def calculate_k_constant_mag(aperture_radius: float, seeing: float) -> float:
    """k(r) once converted into a magnitude."""
    k_constant = calculate_k_constant(aperture_radius, seeing)
    return -2.5 * np.log10(k_constant)
