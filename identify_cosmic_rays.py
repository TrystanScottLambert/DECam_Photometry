"""
Detecting cosmic rays. Going through the individual frames from
the DECam run and seeing if any of the candidates were hit by cosmic rays.
"""

from astropy.io import fits
from astropy.wcs import WCS
import pylab as plt
import numpy as np


infile = '../N964_single_exposures/0050b12beda2b1d1b1306a792ef4cdc4_c4d_210831_073002_ooi_N964_v1.fits.fz'
candidates_file = 'candidates.txt'

candidate_ra, candidate_dec = np.loadtxt(candidates_file, unpack=True)

hdu = fits.open(infile)
wcs = WCS(hdu[1].header)
