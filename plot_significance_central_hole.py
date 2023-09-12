"""
Working out the empirical significance of the central hole according to eduardo.
"""

from rich.progress import track
import numpy as np
from scipy.stats import norm
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import pylab as plt
from plotting import start_plot, end_plot


def count_sources(
        ra_sources: np.ndarray, dec_sources: np.ndarray, x_val: float, y_val: float, radius: float
        ) -> int:
    """
    Counts the number of sources within the search radius at the given center.

    ra_sources, dec_sources = arrays of lyman alpha candidates
    center = searching around given center
    raidus = search radius in degrees.
    """

    candidates = SkyCoord(ra = ra_sources * u.deg, dec = dec_sources * u.deg)
    central_position = SkyCoord(ra = x_val* u.deg,  dec = y_val * u.deg)
    diff = candidates.separation(central_position)
    count = len(np.where(diff < radius*u.deg)[0])
    return count


def mask_central_region(ra_array: np.ndarray, dec_array: np.ndarray, radius: float):
    """
    Removes central apertures that are overlappting with the central area.
    """
    RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
    DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1
    inner_region = 0.271 #deg
    limit = inner_region + radius
    c = SkyCoord(ra_array * u.deg, dec_array * u.deg)
    quasar = SkyCoord(RA_QSO * u.deg, DEC_QSO * u.deg)
    distances = c.separation(quasar)
    cut = np.where(distances > limit * u.deg)
    return cut

def is_valid(mask, center_x: int, center_y: int, radius: int):
    """
    decides whether or not a aperture position is valid based on if the aperture 
    exceeds the limits of the mask.

    center_x, center_y, and radius should be in pixel units.
    """
    cutout = Cutout2D(mask[0].data, (center_x, center_y), radius)
    if np.min(cutout.data) == 0:
        result = False
    elif np.max(cutout.data) ==1:
        result = True
    else:
        print('wtf')
    return result




if __name__ == '__main__':
    infile = 'candidates_e.txt'
    mask_infile = fits.open('DECAM_MASK.fits')
    shape = mask_infile[0].data.shape
    wcs = WCS(mask_infile[0].header)
    pix_scale = mask_infile[0].header['CD2_2']
    ra, dec = np.loadtxt(infile, unpack=True)

    radii = np.arange(0.27, 0.4, 0.02)
    for radius in radii:
        print(radius)
        #radius = 0.27#159338
        
        radius_pixel = radius / pix_scale

        TRIALS = 10000

        center_xs, center_ys = np.random.randint(0, shape[1], TRIALS), np.random.randint(0, shape[0], TRIALS)
        center_ras, center_decs = wcs.pixel_to_world_values(center_xs, center_ys)

        mask = mask_central_region(center_ras, center_decs, radius)
        #center_xs, center_ys, center_ras, center_decs = center_xs[mask], center_ys[mask], center_ras[mask], center_decs[mask]
        valid = np.array([is_valid(mask_infile, x, y, radius_pixel) for x, y in track(zip(center_xs, center_ys))])
        mask = np.where(valid == True)
        numbers = np.array([count_sources(ra, dec, x, y, radius) for x, y in zip(center_ras, center_decs)])

        final_numbers = numbers[valid]

        RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
        DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1
        counts_around_qso = count_sources(ra, dec, RA_QSO, DEC_QSO, radius)

        REDSHIFT_QSO = 6.9
        COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
        ARCSEC_PER_KPC = COSMO.arcsec_per_kpc_proper(REDSHIFT_QSO)
        DEG_PER_MPC = ARCSEC_PER_KPC.to(u.deg / u.Mpc)

        radius_mpc = radius *u.deg / DEG_PER_MPC

        fig = start_plot('# Galaxies within radius', 'Probability')
        ax = fig.add_subplot(111)

        bins = np.arange(0, 11, 1)
        counts, _ = np.histogram(final_numbers, bins = bins)
        x_plotting  = [(bins[i] + bins[i+1])/2 for i in range(len(bins) -1)]

        percents = np.round(counts/len(final_numbers), 2)
        sigmas =  [norm.ppf(prob) for prob in percents]
        # adding second y-axis
        def prob_to_sigma(prob):
            return norm.ppf(prob)
        
        def sigma_to_prob(sigma):
            return norm.cdf(sigma)

        #plt.hist(final_numbers, histtype='step', label=f'radius = {np.round(radius,2)}', bins = np.arange(0, 11, 1))
        ax.step(x_plotting, percents, label = np.round(radius_mpc, 2))
        ax1 = ax.secondary_yaxis('right', functions=(prob_to_sigma, sigma_to_prob))
        ax.axvline(counts_around_qso, ls=':', color='k')
        ax.legend()
        end_plot(f'center_sig_{radius}.png')
    