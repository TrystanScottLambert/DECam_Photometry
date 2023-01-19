""" Plotting the QSO spectra and working out the mag. """

import numpy as np
import pylab as plt


if __name__ == '__main__':
    INFILE = '../QSO_Spectra/J2348_dered.txt'
    wave, flux, err = np.loadtxt(INFILE, unpack=True)
    plt.plot(wave, flux)
    plt.xlabel('Wavelength [angstrom]', fontsize=14)
    plt.ylabel('Flux [Units]', fontsize = 14)
    plt.show()
