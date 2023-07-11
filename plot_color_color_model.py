"""
Modelling color color plot selection.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from synphot import Observation, SpectralElement, SourceSpectrum

from synthetic_lae_spectrum import LaeSpectrum


class Filter:
    """Represantation of a filter."""
    def __init__(self, file_name: str) -> None:
        """Initializing."""
        self.band = SpectralElement.from_file(file_name, wave_unit='nm')

    def observe_magnitude(self, spectrum: SourceSpectrum) -> Observation:
        """Creates an observation through the filter and determines the AB mag."""
        obs = Observation(spectrum, self.band, force='taper')
        return obs.effstim(u.ABmag)


def calculate_color(filter_1: Filter, filter_2: Filter, source: SourceSpectrum) -> float:
    """Calculates the color term for the two given filters of the source."""
    mag_1 = filter_1.observe_magnitude(source)
    mag_2 = filter_2.observe_magnitude(source)
    return mag_1 - mag_2


if __name__ == '__main__':
    Z_BANDPASS_FILE = '../QSO_Spectra/decam_z_bandpass.txt'
    I_BANDPASS_FILE = '../QSO_Spectra/decam_i_bandpass.txt'
    N_BANDPASS_FILE = '../QSO_Spectra/NB964_DECam_F29.txt'
    z_band = Filter(Z_BANDPASS_FILE)
    i_band = Filter(I_BANDPASS_FILE)
    n_band = Filter(N_BANDPASS_FILE)
    plot_redshift_range = np.arange(3, 7.5, 0.01)
    point_redshift_range = np.arange(3, 8, 1)

    def work_out_colors_for_redshift(redshift_range: np.ndarray):
        """Determines the color terms i - z and z-nb for the given redshift range."""
        izs = []
        znbs = []
        for redshift in redshift_range:
            print(redshift)
            spectrum  = LaeSpectrum(redshift).spectrum
            i_z = calculate_color(i_band, z_band, spectrum)
            z_nb = calculate_color(z_band, n_band, spectrum)
            izs.append(i_z)
            znbs.append(z_nb)
        return izs, znbs

    plot_i_color, plot_z_color = work_out_colors_for_redshift(plot_redshift_range)
    point_i_color, point_z_color = work_out_colors_for_redshift(point_redshift_range)
    plt.plot(plot_i_color, plot_z_color)
    plt.scatter(point_i_color, point_z_color, marker='s', color='r')
    plt.show()
