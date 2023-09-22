"""
Basically the same as identify_candidates.py but using custom selection criteria.
"""

from dataclasses import dataclass
import warnings
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
import pylab as plt

import convert_sexcat_into_region as con
from postage_stamps import cut_out_mulitple_stamps
from gui import start_gui
from zero_points import zero_points, ZeroPoints, ZeroPoint
from zero_points_cdfs import zero_points_cdfs
from snr_fit import exponential_func


RANDOM_STATE = 100
np.random.seed(RANDOM_STATE)

# Exponential fit values from snr_fit.py
EXPONENTIAL_FIT_VALS = {
    'i': [90096515422.69316, -0.9205260436876037],
    'z': [84922890028.87346, -0.9206443374807416],
    'n964': [39261983356.719025, -0.9225050055176908],
    'i_cdfs': [218670450252.50992, -0.9078648179670485],
    'z_cdfs': [228699895422.94614, -0.9181554190676867],
    'n964_cdfs': [29312290374.71675, -0.9063544384857863]
}

def calculate_snr(mag_err: float) -> float:
    """Converts the magnitude error into a snr value."""
    return (2.5/np.log(10))/mag_err


def read_all(catalog_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reads in the magnitudes, errors, and works out the snr."""
    mag, err = np.loadtxt(catalog_name, usecols=(4, 5), unpack=True)
    snr = calculate_snr(err)
    return mag, err, snr

def write_region_file(
        ra_array:np.ndarray, dec_array:np.ndarray, outfile:str, size:float = 2.) -> None:
    """Writes a region file with decimal r_a and dec arrays."""
    positions = [con.convert_decimal_degrees_into_celestial(ra_array[i], dec_array[i]) \
                  for i in range(len(ra_array))]

    with open(outfile, 'w', encoding='utf8') as file:
        for pos in positions:
            file.write(f'circle {pos} {size}" # width=4\n')
        file.close()

def write_txt_file(ra_array: np.ndarray, dec_array: np.ndarray, outfile:str) -> None:
    """Writes a text file of decimal r_a dn dec arrays."""
    with open(outfile, 'w', encoding='utf8') as file:
        file.write('# RA DEC \n')
        for i, _ in enumerate(ra_array):
            file.write(f'{ra_array[i]} {dec_array[i]} \n')

def write_output(
        ra_candidates: np.ndarray, dec_candidates: np.ndarray, out_suffix: str, **kwargs) -> None:
    """Writes both a text file and a region file of the given candidates."""
    write_region_file(ra_candidates, dec_candidates, out_suffix+'.reg', **kwargs)
    write_txt_file(ra_candidates, dec_candidates, out_suffix+'.txt')

def update_candidate_red_list(ra_array: np.ndarray, dec_array: np.ndarray, red_list: str) -> None:
    """Updates the red list of candidates which are banned from processing. (obvious artificats)"""
    with open(red_list,'a+', encoding='utf8') as file:
        for i, _ in enumerate(ra_array):
            file.write(f'{ra_array[i]} {dec_array[i]} \n')

def get_red_values(red_list:str) -> tuple[np.ndarray, np.ndarray]:
    """Reads in the r_a and dec of the red sources. (Not allowed to be used)"""
    try:
        ra_rejected, dec_rejected = np.loadtxt(red_list, unpack=True)
    except (OSError, ValueError):
        warnings.warn(f'{red_list} not found or empty. Assuming there are no rejects.')
        ra_rejected, dec_rejected = None, None
    return ra_rejected, dec_rejected

def remove_bad_values(
        ra_array: np.ndarray, dec_array: np.ndarray, red_list:str
        ) -> tuple[np.ndarray, np.ndarray]:
    """Removes the previously rejected r_a and dec values from the current candidates."""
    ra_bad, dec_bad = get_red_values(red_list)
    if ra_bad is not None:
        c_bad = SkyCoord(ra = ra_bad *u.deg, dec = dec_bad * u.deg)
        catalog = SkyCoord(ra = ra_array*u.deg, dec = dec_array*u.deg)
        idx, d2d, _ = c_bad.match_to_catalog_sky(catalog)

        msk = d2d < 1*u.arcsec
        if len(msk) == 1:  # edge case of n=1. Then value isn't read as an array but as a float.
            idx = np.ndarray([idx])
        idx_bad = idx[msk]
        idx_good = np.setdiff1d(np.arange(len(ra_array)), idx_bad)
    else:
        idx_good = np.arange(len(ra_array))


    return ra_array[idx_good], dec_array[idx_good]

@dataclass
class Inputs:
    """Inputs for running the candidate determination."""
    red_list_name: str
    output_name: str
    infile_n964: str
    infile_n964_135: str
    infile_i: str
    infile_z: str
    zero_point_function: ZeroPoints
    images: list[str]
    aperture_radii: float

    @property
    def fits_objects(self):
        """Reads in the images as fits objects."""
        return [fits.open(image) for image in self.images]


class Selection:
    """Selection using the criteria by banados 2013 & Mazzucchelli 2017."""
    def __init__(self, inputs: Inputs, i_2sigma_lim: float, z_2sigma_lim: float):
        self.inputs = inputs
        self.i_lim = i_2sigma_lim
        self.z_lim = z_2sigma_lim

    @property
    def n964_data(self) -> tuple:
        """Reads in the narrowband data and converts instrumental mags into AB mags."""
        inst_mag_n964, inst_mag_n964_err, snr = read_all(self.inputs.infile_n964)
        mag_n964 = inst_mag_n964 +\
              self.inputs.zero_point_function.n964_band.mag_correct(self.inputs.aperture_radii)
        return mag_n964, inst_mag_n964_err, snr

    @property
    def z_data(self) -> tuple:
        """Reads in the z band data and converts instrumental mags into AB mags."""
        inst_mag_z, inst_mag_z_err, z_snr = read_all(self.inputs.infile_z)
        mag_z = inst_mag_z +\
              self.inputs.zero_point_function.z_band.mag_correct(self.inputs.aperture_radii)
        cut = np.where(z_snr < 2)[0]
        mag_z[cut] = self.z_lim
        return mag_z, inst_mag_z_err, z_snr

    @property
    def i_data(self) -> tuple:
        """Reads in the i band data and converts instrumental mags into AB mags."""
        inst_mag_i, inst_mag_i_err, i_snr = read_all(self.inputs.infile_i)
        mag_i = inst_mag_i +\
              self.inputs.zero_point_function.i_band.mag_correct(self.inputs.aperture_radii)
        cut = np.where(i_snr < 2)[0]
        mag_i[cut] = self.i_lim
        return mag_i, inst_mag_i_err, i_snr

    def narrow_color_select(self):
        """
        Looking for excess in narrow band and that that excess is significatnt.
        z-NB> 0.78 and |z-nb| > 2.5 sqrt(u(z)^2 + u(nb)^2)
        """

        z_mag, z_err, _ = self.z_data
        n964_mag, n964_err, _ = self.n964_data
        color = z_mag - n964_mag
        significance = 2.5 * np.hypot(z_err, n964_err)
        first_cut = np.where(color > 0.78)[0]
        final_cut = []
        for idx in first_cut:
            if np.abs(color[idx]) > significance[idx]:
                final_cut.append(idx)
        return np.array(final_cut)

    def continuum_color_select(self):
        """Looking for the continuum break via I-Z > 1.0"""
        z_mag, _, _ = self.z_data
        i_mag, _, _ = self.i_data
        color = i_mag - z_mag
        cut = np.where(color > 1.)[0]
        return cut

    def select_i_band(self):
        """Implementing an I-band cut."""
        _, _, i_snr = self.i_data
        cut = np.where(i_snr < 2)[0]
        return cut

    def apply_selection_criteria(self) -> np.ndarray:
        """Does the selection"""
        narrow_band_excess = self.narrow_color_select()
        continuum_break = self.continuum_color_select()
        i_non_detect = self.select_i_band()
        return np.intersect1d(np.intersect1d(narrow_band_excess, continuum_break), i_non_detect)


class LagerSelection(Selection):
    """
    Class for the Lager selection using the cdfs depth.
    Only important when we make the SNR cut for the i band.
    """

    def __init__(
            self, inputs: Inputs, i_2sigma_lim: float, z_2sigma_lim: float,
            i_2sigma_us: float) -> None:
        super().__init__(inputs, i_2sigma_lim, z_2sigma_lim)
        self.i_lim_us = i_2sigma_us

    def select_i_band(self):
        cut = np.where(self.i_data[0] >= self.i_lim_us)[0]
        return cut

class DegradedLagerSelection(Selection):
    """
    Degrading the CDFS data so that it is equivelent to our depths;
    using the snr fits from snr_fit.py
    """

    def __init__(
            self, inputs: Inputs, i_2sigma_lim: float, z_2sigma_lim: float, n_2sigma_lim: float):
        super().__init__(inputs, i_2sigma_lim, z_2sigma_lim)
        self.n_lim = n_2sigma_lim

    def _prepare_band(
            self, input_file: str, zeropoint: ZeroPoint, sigma_2_depth: float, band: str) -> tuple:
        """Works out the degraded band data."""
        inst_mag, mag_err, snr = read_all(input_file)
        mag = inst_mag + zeropoint.mag_correct(1)
        non_detections = np.where(snr <= 1)[0]
        detections = np.where(snr > 1)[0]
        print('percentage detections', len(non_detections)/len(snr))
        mag[non_detections] = sigma_2_depth
        mag_err[non_detections] = 99.
        old_mag_err_detections = mag_err[detections]
        mag_err[detections] = calculate_snr(
            exponential_func(mag[detections], *EXPONENTIAL_FIT_VALS[band]))
        mag[detections] = np.random.normal(
            mag[detections], np.sqrt(mag_err[detections]**2 - old_mag_err_detections**2))
        snr = calculate_snr(mag_err)
        return mag, mag_err, snr

    @property
    def i_data(self) -> tuple[float, float, float]:
        mag_i, mag_i_err, i_snr = self._prepare_band(
            self.inputs.infile_i, self.inputs.zero_point_function.i_band, self.i_lim, 'i')
        cut = np.where(i_snr < 2)[0]
        mag_i[cut] = self.i_lim
        return mag_i, mag_i_err, i_snr

    @property
    def z_data(self) -> tuple[float, float, float]:
        mag_z, mag_z_err, z_snr = self._prepare_band(
            self.inputs.infile_z, self.inputs.zero_point_function.z_band, self.z_lim, 'z')
        cut = np.where(z_snr < 2)[0]
        mag_z[cut] = self.z_lim
        return mag_z, mag_z_err, z_snr

    @property
    def n964_data(self) -> tuple[float, float, float]:
        mag_n, mag_n_err, n_snr = self._prepare_band(
            self.inputs.infile_n964, self.inputs.zero_point_function.n964_band, self.n_lim, 'n964')
        return mag_n, mag_n_err, n_snr


    def plot_color_color(self) -> None:
        """Plotting the color color plot after degradation and non detecion at our depths"""
        i_z = self.i_data[0] - self.z_data[0]
        z_n = self.z_data[0] - self.n964_data[0]
        plt.scatter(i_z, z_n, s=1)
        plt.axhline(0.78, color='r', ls='--')
        plt.axvline(1, color='r', ls='--')
        plt.xlim(-1, 2.5)
        plt.ylim(-2, 4)
        plt.show()

    def plot_first_color_color(self) -> None:
        """Plotting the color color plot of the pure degradation."""
        mag_i, _, _ = self._prepare_band(
            self.inputs.infile_i, self.inputs.zero_point_function.i_band, self.i_lim, 'i')
        mag_z, _, _ = self._prepare_band(
            self.inputs.infile_z, self.inputs.zero_point_function.z_band, self.z_lim, 'z')
        mag_n, _, _ = self._prepare_band(
            self.inputs.infile_n964, self.inputs.zero_point_function.n964_band, self.n_lim, 'n964')
        i_z = mag_i - mag_z
        z_n = mag_z - mag_n
        plt.scatter(i_z, z_n, s=1)
        plt.axhline(0.78, color='r', ls='--')
        plt.axvline(1, color='r', ls='--')
        plt.xlim(-1, 2.5)
        plt.ylim(-2, 4)
        plt.show()

def perform_selection(selection: Selection):
    """Opens the gui for the user to reject candidates and write to file."""
    cut = selection.apply_selection_criteria()
    print(f'Final candidate count is: {len(cut)}')

    r_a, dec = np.loadtxt(selection.inputs.infile_n964, usecols=(0,1), unpack=True)
    r_a, dec = r_a[cut], dec[cut]
    r_a, dec = remove_bad_values(r_a, dec, selection.inputs.red_list_name)
    print(f'After removing previous rejects, count is: {len(r_a)}')

    i_bands, z_bands, n_bands = cut_out_mulitple_stamps(
        r_a, dec, selection.inputs.fits_objects, pad=60)
    artifacts, candidates = start_gui(i_bands, z_bands, n_bands)
    ra_rejects = r_a[artifacts]
    dec_rejects = dec[artifacts]
    update_candidate_red_list(ra_rejects, dec_rejects, selection.inputs.red_list_name)
    write_output(r_a[candidates], dec[candidates], selection.inputs.output_name, size=8)

if __name__ == '__main__':
    our_inputs = Inputs(
        red_list_name='candidates_red_list.txt',
        output_name='candidates_e',
        infile_n964='../correct_stacks/N964/n964.cat',
        infile_n964_135='../correct_stacks/N964/n964_135.cat',
        infile_i = '../correct_stacks/N964/i.cat',
        infile_z='../correct_stacks/N964/z.cat',
        zero_point_function=zero_points,
        images=(
        '../correct_stacks/N964/i.fits',
        '../correct_stacks/N964/z.fits',
        '../correct_stacks/N964/n964.fits'),
        aperture_radii=1.
    )

    cdfs_inputs = Inputs(
        red_list_name='candidates_red_list_cdfs.txt',
        output_name='candidates_cdfs_e',
        infile_n964='../CDFS_LAGER/n964_cdfs.cat',
        infile_n964_135='../CDFS_LAGER/n964_135_cdfs.cat',
        infile_i = '../CDFS_LAGER/i_cdfs.cat',
        infile_z='../CDFS_LAGER/z_cdfs.cat',
        zero_point_function=zero_points_cdfs,
        images=(
            '../CDFS_LAGER/CDFS_i.fits',
            '../CDFS_LAGER/CDFS_z.fits',
            '../CDFS_LAGER/CDFS_NB.fits'),
        aperture_radii=1.
    )

    I_DEPTH_2_SIGMA = 26.64
    Z_DEPTH_2_SIGMA = 26.58
    N_DEPTH_2_SIGMA = 25.69

    I_DEPTH_2_SIGMA_CDFS = 28.10
    Z_DEPTH_2_SIGMA_CDFS = 27.73

    our_selection = Selection(our_inputs, I_DEPTH_2_SIGMA, Z_DEPTH_2_SIGMA)
    #cdfs_selection = LagerSelection(
    #    cdfs_inputs, I_DEPTH_2_SIGMA_CDFS, Z_DEPTH_2_SIGMA_CDFS, I_DEPTH_2_SIGMA)
    cdfs_selection = DegradedLagerSelection(
        cdfs_inputs, I_DEPTH_2_SIGMA, Z_DEPTH_2_SIGMA, N_DEPTH_2_SIGMA)

    #cdfs_selection.plot_first_color_color()
    #cdfs_selection.plot_color_color()

    #perform_selection(our_selection)
    perform_selection(cdfs_selection)

    #true_cdfs_inputs = cdfs_inputs
    #true_cdfs_inputs.output_name = 'candidates_true_cdfs'
    #true_cdfs_selection = Selection(true_cdfs_inputs, I_DEPTH_2_SIGMA_CDFS, Z_DEPTH_2_SIGMA_CDFS)
    #perform_selection(true_cdfs_selection)
