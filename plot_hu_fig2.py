"""
Script to make figure 2 from Hu et. al., 2019
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from plotting import prettify_axis, end_plot, apply_global_settings
from synthetic_lae_spectrum import LaeSpectrum
from plot_color_color_model import Filter, calculate_color

apply_global_settings()

if __name__ == '__main__':
    REDSHIFT_STEP = 0.001
    redshift_range = np.arange(6.5, 7+REDSHIFT_STEP, REDSHIFT_STEP)
    equivalent_widths = np.array([5, 10, 20, 30, 40, 50])
    Z_BANDPASS_FILE = '../QSO_Spectra/decam_z_bandpass.txt'
    N_BANDPASS_FILE = '../QSO_Spectra/NB964_DECam_F29.txt'
    z_band = Filter(Z_BANDPASS_FILE)
    n_band = Filter(N_BANDPASS_FILE)

    filter_plot_x = n_band.band.to_spectrum1d().spectral_axis.value
    filter_plot_redshift = (filter_plot_x/LaeSpectrum.LLYA_REST.value) - 1
    filter_plot_y = n_band.band.to_spectrum1d().flux.value

    y_values = []
    for ew in equivalent_widths:
        z_nbs = []
        LaeSpectrum.EW_LYA = ew * u.angstrom
        for redshift in redshift_range:
            spectrum = LaeSpectrum(redshift).spectrum
            z_nb = calculate_color(z_band, n_band, spectrum)
            z_nbs.append(z_nb.value)
        y_values.append(np.array(z_nbs))

    fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
    ax_main = fig.add_subplot(111)
    for i, ew in enumerate(equivalent_widths):
        ax_main.plot(redshift_range, y_values[i], label=f'EW = {ew}' + r'$\AA$')
        ax_main.legend(fontsize=12, frameon=False)
        ax_main.set_ylim(0, 3)
        ax_main.set_xlim(6.5, 7.05)
        ax_main.axhline(0.75, ls='--', color='k', alpha=0.5, lw=1.5)
        prettify_axis(ax_main, 'Redshift', 'Z - NB964')

    ax_filter = ax_main.twinx()
    ax_filter.plot(filter_plot_redshift, filter_plot_y, ls=':', color='k', alpha=0.3)
    ax_filter.plot([ax_main.get_xlim()[0],filter_plot_redshift[0]], [0, 0], ls=':', color='k', alpha=0.3)
    ax_filter.fill_between(x= filter_plot_redshift, y1= filter_plot_y,color= "k",alpha= 0.2)
    prettify_axis(ax_filter, '', 'Transmission (%)')
    end_plot('plots/hu_fig2.png')
