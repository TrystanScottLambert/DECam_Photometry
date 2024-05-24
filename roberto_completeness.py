"""
Script to perform the completeness correction as suggested by 
Roberto, in order to address the referee comment.
"""

import numpy as np
import pylab as plt
import astropy.units as u

from synthetic_lae_spectrum import LaeSpectrum
from plot_color_color_model import Filter


"""
Minimum EW based on the color criteria Z-NB > 1
"""



# Reading in the filters.
Z_BANDPASS_FILE = '../QSO_Spectra/decam_z_bandpass.txt'
N_BANDPASS_FILE = '../QSO_Spectra/NB964_DECam_F29.txt'
z_band = Filter(Z_BANDPASS_FILE)
n_band = Filter(N_BANDPASS_FILE)

colors = []
equivalent_widths = np.arange(0.001, 100) #angstroms
for width in equivalent_widths:
	LaeSpectrum.EW_LYA = width * u.angstrom
	LaeSpectrum.FWHM_LYA = 1 * (u.km/u.s)
	lae = LaeSpectrum(6.9)
	color = z_band.observe_magnitude(lae.spectrum) - n_band.observe_magnitude(lae.spectrum)
	colors.append(color.value)

plt.plot(equivalent_widths, colors)
plt.xlabel('EWs', fontsize=30)
plt.ylabel('Z-NB', fontsize=30)
plt.show()


# Maxmimum EW based on the depth plots and z-detection criteria
# |z-NB| > 2.5*sqrt(sigma_z**2 + sigma_NB**2)

