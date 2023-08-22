"""
Script to work out the statistical signficance of the inner region. 
"""

from rich.progress import track
import numpy as np
import pylab as plt
from scipy.stats import poisson, norm
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
import astropy.units as u
from plot_radial_distribution import Mask


REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
ARCSEC_PER_KPC = COSMO.arcsec_per_kpc_comoving(REDSHIFT_QSO)
DEG_PER_MPC = ARCSEC_PER_KPC.to(u.deg / u.Mpc)
RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24) * u.deg
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1 *u.deg
qso_position = SkyCoord(ra = RA_QSO, dec = DEC_QSO)

def calculate_distances_to_quasar(ra_array: np.ndarray, dec_array: np.ndarray):
    """Works out the distances to the quasar for a given ra and dec range."""
    coords = SkyCoord(ra= ra_array*u.deg, dec = dec_array*u.deg)
    return qso_position.separation(coords)

def count_inner_outer(ra_array: np.ndarray, dec_array: np.ndarray, distance_deg: float) -> tuple:
    """Calculates the number of sources within the deg distance and how many are outside of it."""
    distances = calculate_distances_to_quasar(ra_array, dec_array)
    results = distances < distance_deg
    inner = len(np.where(results == True)[0])
    outer = len(np.where(results == False)[0])

    return inner, outer

def calculate_area_ratio(mask: Mask, radius: float) -> float:
    """Works out the scaling factor for the counts."""
    radius_pixels = radius / mask.deg_per_pix
    center = mask.wcs.world_to_pixel_values(RA_QSO.value, DEC_QSO.value)
    outer_area = mask.calculate_area(center, radius_pixels, 100000)
    inner_area = mask.calculate_area(center, 0, radius_pixels)
    ratio = inner_area/outer_area
    return ratio.value

if __name__ == '__main__':
    decam = Mask('DECAM_MASK.fits')
    ra, dec = np.loadtxt('candidates_e.txt', unpack=True)

    distances = np.arange(70, 76, 1) *u.Mpc

    sigmas = []
    for distance in track(distances):
        inner_region_distance_deg = DEG_PER_MPC * distance
        inner, outer = count_inner_outer(ra, dec, inner_region_distance_deg)
        ratio = calculate_area_ratio(decam, inner_region_distance_deg)
        expected_inner = outer * ratio

        probability = poisson.cdf(inner, expected_inner)
        sigmas.append(norm.ppf(probability))

    plt.plot(distances, sigmas)
    plt.show()
