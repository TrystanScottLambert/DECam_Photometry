"""
Basically the same as identify_candidates.py but using custom selection criteria.
"""

from dataclasses import dataclass
from abc import abstractmethod
from typing import Protocol
import warnings
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits

import convert_sexcat_into_region as con
import postage_stamps as ps
from gui import start_gui
from zero_points import zero_points, ZeroPoints
from zero_points_cdfs import zero_points_cdfs


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
    """Writes a region file with decimal ra and dec arrays."""
    positions = [con.convert_decimal_degrees_into_celestial(ra_array[i], dec_array[i]) \
                  for i in range(len(ra_array))]
    file = open(outfile, 'w', encoding='utf8')
    for pos in positions:
        file.write(f'circle {pos} {size}" # width=4\n')
    file.close()

def write_txt_file(ra_array: np.ndarray, dec_array: np.ndarray, outfile:str) -> None:
    """Writes a text file of decimal ra dn dec arrays."""
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
    """Reads in the ra and dec of the red sources. (Not allowed to be used)"""
    try:
        ra_rejected, dec_rejected = np.loadtxt(red_list, unpack=True)
        print(ra_rejected)
    except (OSError, ValueError):
        warnings.warn(f'{red_list} not found or empty. Assuming there are no rejects.')
        ra_rejected, dec_rejected = None, None
    return ra_rejected, dec_rejected

def remove_bad_values(
        ra_array: np.ndarray, dec_array: np.ndarray, red_list:str
        ) -> tuple[np.ndarray, np.ndarray]:
    """Removes the previously rejected ra and dec values from the current candidates."""
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


class SelectionCriteria(Protocol):
    """Abstract base class for selection criteria."""

    inputs: Inputs
    n964_2_lim: float
    n964_135_lim: float
    i_lim: float
    z_lim: float

    @abstractmethod
    def select_n964(self) -> np.ndarray:
        """The selection criteria for the n964 filter."""

    @abstractmethod
    def select_i(self) -> np.ndarray:
        """The selection criteria for the i band filter."""

    @abstractmethod
    def select_z(self) -> np.ndarray:
        """The selection criteria for the z band filter."""

    def apply_selection_criteria(self) -> np.ndarray:
        """Applies all the selction criteria to the input."""
        msk_n964 = self.select_n964()
        msk_i = self.select_i()
        msk_z = self.select_z()
        first = np.intersect1d(msk_n964, msk_i)
        return np.intersect1d(first, msk_z)

    @property
    @abstractmethod
    def n964_data(self) -> tuple:
        """Read in the n964 data."""

    @property
    @abstractmethod
    def i_data(self) -> tuple:
        """Read in the n964 data."""

    @property
    @abstractmethod
    def z_data(self) -> tuple:
        """Read in the n964 data."""


class MagCutSelection(SelectionCriteria):
    """Selection using mag cuts."""

    def __init__(
            self, inputs: Inputs, n964_2_lim: float, n964_135_lim: float, i_lim: float, z_lim: float
            ) -> None:
        self.inputs = inputs
        self.n964_2_lim = n964_2_lim
        self.n964_135_lim =  n964_135_lim
        self.i_lim = i_lim
        self.z_lim = z_lim

    @property
    def n964_data(self) -> tuple:
        inst_mag_n964, _, _ = read_all(self.inputs.infile_n964)
        inst_mag_135, _, _ = read_all(self.inputs.infile_n964_135)
        mag_n964 = inst_mag_n964 + self.inputs.zero_point_function.n964_band.mag_correct(
            self.inputs.aperture_radii)
        mag_n964_135 = inst_mag_135 + self.inputs.zero_point_function.n964_band.mag_correct(
            1.35/2)
        return mag_n964, mag_n964_135

    def select_n964(self) -> np.ndarray:
        #1. mag < 24.2 for the N964 filter in 2" and mag < 24 in 1.35" apertures.
        two_arcsecond_cut = np.where(self.n964_data[0] < self.n964_2_lim)[0]
        one_arcsecond_cut = np.where(self.n964_data[1] < self.n964_135_lim)[0]
        n964_cut = np.intersect1d(two_arcsecond_cut, one_arcsecond_cut)
        return n964_cut

    @property
    def i_data(self) -> tuple:
        inst_mag_i, _, _ = read_all(self.inputs.infile_i)
        mag_i = inst_mag_i + self.inputs.zero_point_function.i_band.mag_correct(
            self.inputs.aperture_radii)
        return (mag_i,)

    def select_i(self) -> np.ndarray:
        #2. mag > 25.8 for the i band filter.
        return np.where(self.i_data[0] > self.i_lim)[0]

    @property
    def z_data(self) -> tuple:
        inst_mag_z, _, _ = read_all(self.inputs.infile_z)
        mag_z = inst_mag_z + self.inputs.zero_point_function.z_band.mag_correct(
            self.inputs.aperture_radii)
        return (mag_z,)

    def select_z(self) -> np.ndarray:
        #3a.  z-N964 > 1.9 and mag < 25.6 for the z filter
        #3.b  S/N_2" > 25.6 for the z filter.
        color_selection = 1.9
        color = self.z_data[0] - self.n964_data[0]

        z_cut_a_1 = np.where(color > color_selection)
        z_cut_a_2 = np.where(self.z_data[0] < self.z_lim)[0]
        z_cut_a = np.intersect1d(z_cut_a_1, z_cut_a_2)

        z_cut_b = np.where(self.z_data[0] > self.z_lim)[0]
        return np.union1d(z_cut_a, z_cut_b)

def perform_selection(selection: MagCutSelection):
    """Opens the gui for the user to reject candidates and write to file."""
    cut = selection.apply_selection_criteria()
    print(cut)
    print(f'Final candidate count is: {len(cut)}')

    ra, dec = np.loadtxt(selection.inputs.infile_n964, usecols=(0,1), unpack=True)
    print(ra)
    ra, dec = ra[cut], dec[cut]
    ra, dec = remove_bad_values(ra, dec, selection.inputs.red_list_name)
    print(f'After removing previous rejects, count is: {len(ra)}')

    i_bands, z_bands, n_bands = [], [], []
    for i, _ in enumerate(ra): 
        i_filter, z_filter, n964_filter = ps.cut_out_stamps(
            ra[i], dec[i], selection.inputs.fits_objects, pad=60)
        i_bands.append(i_filter)
        z_bands.append(z_filter)
        n_bands.append(n964_filter)

    artifacts, candidates = start_gui(i_bands, z_bands, n_bands)
    ra_rejects = ra[artifacts]
    dec_rejects = dec[artifacts]
    update_candidate_red_list(ra_rejects, dec_rejects, selection.inputs.red_list_name)
    write_output(ra[candidates], dec[candidates], selection.inputs.output_name, size=8)

if __name__ == '__main__':
    our_inputs = Inputs(
        red_list_name='candidates_red_list.txt',
        output_name='candidates',
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
        output_name='candidates_cdfs',
        infile_n964='../CDFS_LAGER/n964_cdfs.cat',
        infile_n964_135='../CDFS_LAGER/n964_135_cdfs.cat',
        infile_i = '../CDFS_LAGER/i_cdfs.cat',
        infile_z='../CDFS_LAGER/z_cdfs.cat',
        zero_point_function=zero_points_cdfs,
        images=(
            '../CDFS_LAGER/i.fits',
            '../CDFS_LAGER/z.fits',
            '../CDFS_LAGER/n964.fits'),
        aperture_radii=1.
    )

    our_selection = MagCutSelection(our_inputs, 24.2, 24, 25.8, 25.6)
    cdfs_selection = MagCutSelection(cdfs_inputs, 24.2, 24, 25.8, 25.6)

    perform_selection(our_selection)
    perform_selection(cdfs_selection)
