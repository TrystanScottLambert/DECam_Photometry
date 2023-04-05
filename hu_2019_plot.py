"""Python script to display the DECAM filters."""

from dataclasses import dataclass
import pylab as plt
import numpy as np
import matplotlib as mpl
from astropy.coordinates import SkyCoord
import astropy.units as u
import plotting



mpl.rcParams.update({'font.size': 2})
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10)
mpl.rc('ytick', labelsize=10)


def read_in_cats(filt: str):
    """read in the test cats for DECCAM"""
    r_a, dec, mag, mag_err = np.loadtxt(FILTERS[filt].file_name, unpack=True, usecols=(0, 1, 4, 5))
    non_detections = np.where(mag == 99)[0]
    detections = np.where(mag != 99)[0]

    mag[detections] = mag[detections] + FILTERS[filt].zpt
    mag_err[non_detections] = (2.5/np.log(10)) / 2 #factor of two is just for aesthitics

    return r_a, dec , mag, mag_err


@dataclass
class Filter:
    """Stores the filter information."""
    zpt: float
    limiting_mag: float
    file_name: str
    hu_snr_limit: float


@dataclass
class Quasar:
    """All the magnitude information"""
    n964_mag: float
    z_mag: float
    i_mag: float

    def color(self, filt: str):
        """Determines the color depending on the filter."""
        if filt == 'i':
            color = self.i_mag - self.n964_mag
        elif filt == 'z':
            color = self.z_mag - self.n964_mag
        else:
            raise ValueError('Filter must be either "z" or "i".')
        return color


@dataclass
class ErrorPlotValues:
    """x,y,xerr,yerr values for a errobar plot"""
    x_data: np.array
    y_data: np.array
    x_err: np.array
    y_err: np.array


@dataclass
class MagValues:
    """Storing all the information for the hu plot"""
    narrow_mag: np.array
    broad_mag: np.array
    narrow_err: np.array
    broad_err: np.array
    idx_gals: np.array
    non_detections: np.array

    @property
    def color(self):
        """Calculates the color"""
        return self.broad_mag - self.narrow_mag

    @property
    def color_err(self):
        """Calculates the error of the color"""
        return np.hypot(self.narrow_err, self.broad_err)

    @property
    def gals(self):
        """Returns a ErrorPlotValues object for all the candidates"""
        return ErrorPlotValues(self.narrow_mag[self.idx_gals], self.color[self.idx_gals],
                               self.narrow_err[self.idx_gals], self.color_err[self.idx_gals])

    @property
    def non_detects(self):
        """Determines the values for a errorbar plot of the non detectionss"""
        return ErrorPlotValues(self.narrow_mag[self.idx_gals][self.non_detections],
                               self.color[self.idx_gals][self.non_detections],
                               self.narrow_err[self.idx_gals][self.non_detections],
                               self.color_err[self.idx_gals][self.non_detections])


@dataclass
class PlotSettings:
    """settings for the plots"""
    x_label: str
    y_label: str
    outfile: str


def plot_hu_plot(conf:MagValues, qso:Quasar, plot_conf: PlotSettings, filt:str):
    """makes the hu plot specifically"""
    plotting.start_plot(plot_conf.x_label, plot_conf.y_label)
    plt.scatter(conf.narrow_mag, conf.color, s=1, color='k', alpha=0.5)
    plt.scatter(qso.n964_mag, qso.color(filt), marker='*', s=50, color='m', label='QSO')

    plt.errorbar(conf.gals.x_data, conf.gals.y_data, color='r',
                xerr=conf.gals.x_err, yerr=conf.gals.y_err,
                fmt='o', ms=2, elinewidth=1, label = 'Candidates')

    plt.errorbar(conf.non_detects.x_data,
                FILTERS[filt].limiting_mag - conf.non_detects.x_data, color='r',
                xerr=conf.non_detects.x_err, yerr=conf.non_detects.y_err, fmt='o', ms=2,
                elinewidth=1, lolims=True)

    plt.xlim(19, 25)
    plt.ylim(-1.1,6)
    plt.axhline(FILTERS[filt].hu_snr_limit, color='r', lw=1, ls='--')
    plt.legend(fontsize=8, loc=2)
    plotting.end_plot(plot_conf.outfile)


FILTERS = {
    'z': Filter(30.538, 25.5, '../correct_stacks/N964/z.cat', 1.9),
    'n964': Filter(28.97, 25.0, '../correct_stacks/N964/n964.cat', 0),
    'i': Filter(30.870, 25.5, '../correct_stacks/N964/i.cat', 0.8)
}

RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1

if __name__ == '__main__':
    CANDIDATE_FILE = 'candidates.txt'

    i_ra, i_dec, i_mag, i_mag_err = read_in_cats('i')
    z_ra, z_dec, z_mag, z_mag_err = read_in_cats('z')
    n_ra, n_dec, n_mag, n_mag_err = read_in_cats('n964')

    catalog = SkyCoord(ra = n_ra*u.deg, dec = n_dec*u.deg)
    c_qso = SkyCoord(ra = RA_QSO*u.deg, dec = DEC_QSO*u.deg)
    idx, _,  _ = c_qso.match_to_catalog_sky(catalog)
    our_qso = Quasar(n_mag[idx], z_mag[idx], i_mag[idx])

    gals_ra, gals_dec = np.loadtxt(CANDIDATE_FILE, unpack=True)
    c_gals = SkyCoord(ra = gals_ra*u.deg, dec = gals_dec*u.deg)
    idx_gals, d2d_gals, _ = c_gals.match_to_catalog_sky(catalog)

    non_detections_z = np.where(z_mag[idx_gals] == 99)[0]
    z_mag_settings = MagValues(n_mag, z_mag, n_mag_err, z_mag_err, idx_gals, non_detections_z)
    z_plot_settings = PlotSettings('N964 [Mag]',  'z - N964 [Mag]', 'plots/hu_plot_z.png')
    plot_hu_plot(z_mag_settings, our_qso, z_plot_settings, 'z')

    non_detections_i = np.where(i_mag[idx_gals] == 99)[0]
    i_mag_settings = MagValues(n_mag, i_mag, n_mag_err, i_mag_err, idx_gals, non_detections_i)
    i_plot_settings = PlotSettings('N964 [Mag]',  'i - N964 [Mag]', 'plots/hu_plot_i.png')
    plot_hu_plot(i_mag_settings, our_qso, i_plot_settings, 'i')
