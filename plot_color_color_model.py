"""
Modelling color color plot selection.
"""

from rich.progress import track
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from synphot import Observation, SpectralElement, SourceSpectrum, Empirical1D, etau_madau
import synphot.units as su

from plotting import start_plot, end_plot, prettify_axis
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
    plot_redshift_range = np.arange(5, 8, 0.01)
    point_redshift_range = np.arange(5, 9, 1)

    def work_out_colors_for_redshift(redshift_range: np.ndarray, offset_velocity: float = 0):
        """Determines the color terms i - z and z-nb for the given redshift range."""
        izs = []
        znbs = []
        for redshift in track(redshift_range):
            offset_redshift = redshift - offset_velocity/3e5
            spectrum  = LaeSpectrum(offset_redshift).spectrum
            i_z = calculate_color(i_band, z_band, spectrum)
            z_nb = calculate_color(z_band, n_band, spectrum)
            izs.append(i_z.value)
            znbs.append(z_nb.value)
        return izs, znbs

    VEL_OFFSET = -1000 #km/s
    plot_i_color, plot_z_color = work_out_colors_for_redshift(plot_redshift_range, VEL_OFFSET)
    point_i_color, point_z_color = work_out_colors_for_redshift(point_redshift_range, VEL_OFFSET)
    qso_i_color, qso_z_color = work_out_colors_for_redshift([6.9], VEL_OFFSET)

    start_plot('i - z', 'z - NB964')
    plt.plot(plot_i_color, plot_z_color, lw=2)
    plt.scatter(point_i_color, point_z_color, marker='s', color='r')
    for i_color, z_color, redshift in zip(point_i_color, point_z_color, point_redshift_range):
        plt.text(i_color, z_color, f'{redshift}', fontsize=20, zorder=100)

    plt.axvline(1, ls=':', color='k', alpha=0.5)
    plt.axhline(0.78, ls=':', color='k', alpha=0.5)
    print(qso_i_color, qso_z_color)
    plt.scatter(qso_i_color[0], qso_z_color[0], s=50, marker='*', color='m', zorder=99)
    plt.vlines(x=1, ymin=0.78, ymax=6, color='g', lw=2, zorder=99)
    plt.hlines(y=0.78, xmin=1, xmax=6, color='g', lw=2, zorder=99)
    plt.xlim(-4, 6)
    plt.ylim(-2.5, 6)
    end_plot('plots/tracks.png')
    plt.show()

    i_band_spectrum = i_band.band.to_spectrum1d()
    z_band_spectrum = z_band.band.to_spectrum1d()
    n_band_spectrum = n_band.band.to_spectrum1d()

    offset_red = 6.9018 + VEL_OFFSET/3e5
    spectrum = LaeSpectrum(offset_red).spectrum
    spectrum_1d = spectrum.to_spectrum1d()
    fig = plt.figure(figsize = (1.5 * 3.54, 3.54), dpi = 600)
    ax = fig.add_subplot(111)
    ax.plot(i_band_spectrum.spectral_axis.value/10, i_band_spectrum.flux.value, color='b', ls='--',lw=1, label='DECam-i')
    ax.plot(z_band_spectrum.spectral_axis.value/10, z_band_spectrum.flux.value, color='r', ls=':',lw=1, label='DECam-z')
    #ax.fill_between(x=i_band_spectrum.spectral_axis.value/10 , y1=i_band_spectrum.flux.value,color= "b",alpha= 0.1)
    #ax.fill_between(z_band_spectrum.spectral_axis.value/10, z_band_spectrum.flux.value, color='r', ls=':',alpha=0.1)
    
    ax.fill_between(n_band_spectrum.spectral_axis.value/10, n_band_spectrum.flux.value, color='g',alpha=0.1, lw=1, label='NB964')
    ax.set_xlim(670, 1050)
    #ax.axvline( 91.1267*(7.9), color='k', ls=':')
    prettify_axis(ax, 'Wavelength [nm]', 'Transmission')
    ax_real = ax.twinx()
    ax_real.plot(spectrum_1d.spectral_axis.value/10,spectrum_1d.flux.value, color='k',alpha=0.8, lw=1, label='LAE spectrum')
    ax.plot(n_band_spectrum.spectral_axis.value/10, n_band_spectrum.flux.value, color='g',lw=1)
    ax_real.set_yticks([])
    ax_real.set_yticklabels([])
    #ax_real.legend()
    ax.legend(loc=2)
    end_plot('plots/Filters.png')
    plt.show()
    

    #GIF Chiara
    '''i_band_spectrum = i_band.band.to_spectrum1d()
    z_band_spectrum = z_band.band.to_spectrum1d()
    n_band_spectrum = n_band.band.to_spectrum1d()

    i_colors = []
    z_colors = []
    for i, redshift in enumerate(plot_redshift_range):
        spectrum = LaeSpectrum(redshift).spectrum
        spectrum_1d = spectrum.to_spectrum1d()
        fig = plt.figure()
        ax_track = fig.add_subplot(121)
        ax_spec = fig.add_subplot(122)
        i_mag = i_band.observe_magnitude(spectrum)
        z_mag = z_band.observe_magnitude(spectrum)
        n_mag = n_band.observe_magnitude(spectrum)
        i_colors.append(i_mag.value  - z_mag.value)
        z_colors.append(z_mag.value  - n_mag.value)
        ax_spec.plot(i_band_spectrum.spectral_axis.value, i_band_spectrum.flux.value)
        ax_spec.plot(z_band_spectrum.spectral_axis.value, z_band_spectrum.flux.value)
        ax_spec.plot(n_band_spectrum.spectral_axis.value, n_band_spectrum.flux.value)
        ax_spec.set_xlabel('Wavelength [angstrom]')
        ax_spec.set_ylabel('Filter Transmission')
        ax_spec.set_xlim(5000, 12000)
        ax_track.plot(i_colors, z_colors)
        ax_track.set_xlim(-2, 6)
        ax_track.set_ylim(-3, 3)
        ax_track.axhline(0.75, ls=':', color='k')
        ax_track.axvline(1, ls=':', color='k')
        ax_spec_real = ax_spec.twinx()
        ax_spec_real.plot(spectrum_1d.spectral_axis.value,
                          spectrum_1d.flux.value, color='k',
                          label=f'z = {round(redshift, 1)}')
        ax_spec_real.legend(handlelength=0, handletextpad=0, frameon=False)
        ax_spec_real.set_ylabel('Flux [PHOTLAM]')
        ax_track.set_ylabel('Z - NB965')
        ax_track.set_xlabel('I - Z')
        plt.savefig(f'test_delete_{i}.png')
        plt.close()'''


    #GIF Shapely
    '''infile_shapely = 'shapely_template.dat'
    wave, flux = np.loadtxt(infile_shapely, unpack=True)
    wave = wave*u.angstrom
    flux = flux*su.FNU
    shapely_spectrum = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux, keep_neg=True)

    i_colors = []
    z_colors = []
    for i, redshift in enumerate(np.arange(4, 8, 0.05)):
        spectrum = SourceSpectrum(shapely_spectrum.model, z=redshift, keep_neg=True)
        extinction_curve = etau_madau(spectrum.to_spectrum1d().spectral_axis, redshift)
        spectrum = spectrum * extinction_curve
        spectrum_1d = spectrum.to_spectrum1d()
        fig = plt.figure()
        ax_track = fig.add_subplot(121)
        ax_spec = fig.add_subplot(122)
        i_mag = i_band.observe_magnitude(spectrum)
        z_mag = z_band.observe_magnitude(spectrum)
        n_mag = n_band.observe_magnitude(spectrum)
        i_colors.append(i_mag.value  - z_mag.value)
        z_colors.append(z_mag.value  - n_mag.value)
        ax_spec.plot(i_band_spectrum.spectral_axis.value, i_band_spectrum.flux.value)
        ax_spec.plot(z_band_spectrum.spectral_axis.value, z_band_spectrum.flux.value)
        ax_spec.plot(n_band_spectrum.spectral_axis.value, n_band_spectrum.flux.value)
        ax_spec.set_xlabel('Wavelength [angstrom]')
        ax_spec.set_ylabel('Filter Transmission')
        ax_spec.set_xlim(5000, 12000)
        ax_track.plot(i_colors, z_colors)
        ax_track.set_xlim(-2, 6)
        ax_track.set_ylim(-2.6, 3)
        ax_track.axhline(0.75, ls=':', color='k')
        ax_track.axvline(1, ls=':', color='k')
        ax_spec_real = ax_spec.twinx()
        ax_spec_real.plot(spectrum_1d.spectral_axis.value,
                          spectrum_1d.flux.value, color='k',
                          label=f'z = {round(redshift, 1)}')
        ax_spec_real.legend(handlelength=0, handletextpad=0, frameon=False)
        ax_spec_real.set_ylabel('Flux [PHOTLAM]')
        ax_track.set_ylabel('Z - NB965')
        ax_track.set_xlabel('I - Z')
        plt.savefig(f'test_delete_{i}.png')
        plt.close()'''
