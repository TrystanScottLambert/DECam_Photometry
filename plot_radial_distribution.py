"""
Plotting the radial distribution of the counts.
"""

from rich.progress import track
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import Angle
import pylab as plt
import plotting

RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1
REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
ARCSEC_PER_KPC = COSMO.arcsec_per_kpc_proper(REDSHIFT_QSO)
DEG_PER_MPC = ARCSEC_PER_KPC.to(u.deg / u.Mpc)

class Mask:
    """Mask class"""
    def __init__(self, file_name: str) -> None:
        self.hdu = fits.open(file_name)
        self.data = self.hdu[0].data
        self.wcs = WCS(self.hdu[0].header)
        self.deg_per_pix = self.hdu[0].header['PC2_2']*u.deg
        self.pix_area_deg = (self.deg_per_pix) ** 2

    def calculate_area(self, center: tuple, inner_radius: float, outer_radius: float) -> float:
        """
        Calculates the area by counting the number of non-masked pixels and multiplying
        by the pixel area.
        """
        cutout = create_annulus_cutout(self.data, center, inner_radius, outer_radius)
        number_pixels = len(np.where(cutout==1)[0])
        return number_pixels * self.pix_area_deg

    @property
    def center(self):
        """Calculates the center RA/DEC of the center of the image."""
        rows, columns = self.data.shape
        center_row = rows/2
        center_column = columns/2 
        return self.wcs.pixel_to_world_values(center_column, center_row)

def create_annulus_cutout(image: np.ndarray, center: tuple, inner_radius: float, outer_radius: float):
    """
    Create a new image with an annulus cutout from the input image.
    
    Parameters:
        image (numpy.ndarray): 2D fits image array.
        center (tuple): Center coordinates (y, x) of the annulus.
        inner_radius (float): Inner radius of the annulus.
        outer_radius (float): Outer radius of the annulus.
        
    Returns:
        numpy.ndarray: New image with pixels outside the annulus set to NaN.
    """
    y_indicies, x_indicies = np.indices(image.shape)

    distance_from_center = np.sqrt((x_indicies - center[1])**2 + (y_indicies - center[0])**2)

    mask = np.logical_and(distance_from_center >= inner_radius, distance_from_center <= outer_radius)

    cutout_image = image.copy()
    cutout_image[~mask] = np.nan
    return cutout_image

def calculate_distances_of_candidates(center: tuple, candidate_list: str) -> Angle:
    """
    Works out the distance to the candidates from the given center.
    Center must be in decimal ra and dec and candidate list is the file name
    with the positions of the candidates. 
    """
    center_position = SkyCoord(ra = center[0]*u.deg, dec = center[1]*u.deg)
    ras, decs = np.loadtxt(candidate_list, unpack=True)
    candidates = SkyCoord(ra = ras*u.deg, dec = decs*u.deg)
    separations = center_position.separation(candidates)
    return separations

def plot_radial_profile(counts: np.ndarray, areas_vals: list, distances: np.ndarray, separations: np.ndarray, outfile: str):
    """plots the radial profile"""
    areas = np.array([area.value for area in areas_vals])
    null_count_values = np.where(counts==0)[0]
    non_null_count_values = np.where(counts!=0)[0]
    counts[null_count_values] = 0#2 # poisson distribution means a zero count is less than 2 counts.
    y = counts/areas
    y_err = np.sqrt(counts)/areas


    fig = plt.figure(figsize = (3.54, 3.54/2), dpi = 600)
    ax = fig.add_subplot(111)
    ax.errorbar(distances[non_null_count_values], y[non_null_count_values], yerr = y_err[non_null_count_values], fmt='ok', ecolor='r', ms = 2, capsize=2)
    ax.errorbar(distances[null_count_values], y[null_count_values], yerr = y_err[null_count_values], fmt='ok', ecolor='r', ms = 2, capsize=2, uplims=True)
    ax.set_xlabel('Distance from QSO [pMpc]')
    ax.set_ylabel(r'Surface Density [deg^{-2}$]')
    ax.minorticks_on()
    ax.tick_params(which='both', width=1.2,direction='in')
    ax.tick_params(which='major', length=3, direction='in')
    #ax.set_yscale('log')

    ax1 = ax.twiny()
    ax1.errorbar(separations, y, yerr=y_err, fmt='ok', alpha=0)
    ax1.set_xlabel('Distance from QSO [deg]')
    ax1.minorticks_on()
    ax1.tick_params(which='both', width=1.2,direction='in')
    ax1.tick_params(which='major', length=3, direction='in')
    plotting.end_plot(outfile)
    

if __name__ == '__main__':
    CDFS_MASK_FILE = 'CDFS_MASK.fits'
    DECAM_MASK_FILE = 'DECAM_MASK.fits'

    CDFS_CANDIDATES = 'candidates_cdfs_e.txt'
    DECAM_CANDIDATES = 'candidates_e.txt'

    cdfs = Mask(CDFS_MASK_FILE)
    cdfs_center = (cdfs.data.shape[1]/2, cdfs.data.shape[0]/2)
    decam = Mask(DECAM_MASK_FILE)
    decam_center = decam.wcs.world_to_pixel_values(RA_QSO, DEC_QSO)


    radii_mpc = np.arange(0,21,1) * u.Mpc
    radii_deg = radii_mpc * DEG_PER_MPC
    radii_pix = radii_deg / cdfs.deg_per_pix
    average_radii_mpc = np.array([(radii_mpc[i].value + radii_mpc[i+1].value)/2 for i in range(len(radii_mpc) -1)]) * u.Mpc
    average_radii_deg = average_radii_mpc * DEG_PER_MPC

    areas_cdfs = [cdfs.calculate_area(cdfs_center, radii_pix[i], radii_pix[i+1]) for  i in track(range(len(radii_pix) -1), 'area cdfs')]
    areas_decam = [decam.calculate_area(decam_center, radii_pix[i], radii_pix[i+1]) for  i in track(range(len(radii_pix) -1), 'area decam')]

    distances_cdfs = calculate_distances_of_candidates(cdfs.center, CDFS_CANDIDATES)
    distances_decam = calculate_distances_of_candidates((RA_QSO, DEC_QSO), DECAM_CANDIDATES)

    counts_cdfs = np.array(
        [len(np.where((distances_cdfs.value > radii_deg[i].value) & (distances_cdfs.value < radii_deg[i+1].value))[0]) \
          for i in range(len(radii_deg) -1)])

    counts_decam = np.array(
        [
            len(np.where((distances_decam.value > radii_deg[i].value) & (distances_decam.value < radii_deg[i+1].value))[0]) \
          for i in range(len(radii_deg) -1)
        ]
    )

    plot_radial_profile(counts_cdfs, areas_cdfs, average_radii_mpc, average_radii_deg, 'plots/surface_density_cdfs.png')
    plot_radial_profile(counts_decam, areas_decam, average_radii_mpc, average_radii_deg, 'plots/surface_density.png')
