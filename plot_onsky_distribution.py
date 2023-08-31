"""Plotting the onsky distribution and the distance from the qso distribution"""

import numpy as np
from matplotlib import pyplot as plt
from regions import PixCoord, CirclePixelRegion
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval

from regionfy_catalog import load_region
import plotting


def get_radius(physical_radius_mpc: float) -> float:
    """Given a physical radius in Mpc, return the radius in arcseconds"""
    physical_radius_deg = physical_radius_mpc * u.Mpc * DEG_PER_MPC
    return physical_radius_deg.value

def set_region(ra_pix, dec_pix, radius_mpc: float):
    """Creates a circular region object with the given radius."""
    region = CirclePixelRegion(center=PixCoord(x=ra_pix, y=dec_pix),
                                radius=get_radius(radius_mpc) / DEG_PER_PIX)
    return region


DECAM_SHAPE_FILE = '../correct_stacks/N964/n964_weight.fits'
CDFS_IMAGE_FILE = '../correct_stacks/N964/n964.fits'

image_data = fits.open(CDFS_IMAGE_FILE)
decam_hdu = fits.open(DECAM_SHAPE_FILE)
decam_wcs = WCS(decam_hdu[0].header)

RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1
REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
ARCSEC_PER_KPC = COSMO.arcsec_per_kpc_comoving(REDSHIFT_QSO)
DEG_PER_MPC = ARCSEC_PER_KPC.to(u.deg / u.Mpc)
DEG_PER_PIX = np.abs(decam_hdu[0].header['PC2_2'])
REGION_FILE = 'decam.reg'
region_decam_fov = load_region(REGION_FILE)
region_decam_mask = fits.open('DECAM_MASK.fits')
region_decam_mask = region_decam_mask[0].data
region_decam_mask[region_decam_mask == 0] = np.nan

image_data[0].data = image_data[0].data * region_decam_mask

if __name__ == '__main__':
    INFILE = 'candidates_e.txt'
    ra, dec = np.loadtxt(INFILE, unpack=True)

    ra_plot, dec_plot = decam_wcs.world_to_pixel_values(ra, dec)
    ra_qso_plot, dec_qso_plot = decam_wcs.world_to_pixel_values(RA_QSO, DEC_QSO)

    #On sky distribution plot.
    #region_1 = set_region(ra_qso_plot, dec_qso_plot, 35)
    #region_10 = set_region(ra_qso_plot, dec_qso_plot, 55.2)
    region_20 = set_region(ra_qso_plot, dec_qso_plot, 75)

    fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
    ax = fig.add_subplot(projection = decam_wcs)
    zscale = ZScaleInterval()
    lower, upper = zscale.get_limits(image_data[0].data)
    ax.imshow(image_data[0].data, cmap='gray_r', vmin=lower, vmax=upper)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.imshow(decam_hdu[0].data, alpha=0)
    #region_1.plot(ax=ax, color='red', lw=2.0, ls=':')
    #region_10.plot(ax=ax, color='red', lw=2.0, ls=':')
    region_20.plot(ax=ax, color='red', lw=2.0, ls=':')
    region_decam_fov.plot(ax = ax, color='k', lw=2.0)
    ax.scatter(ra_plot, dec_plot, s=10)
    ax.scatter(ra_qso_plot, dec_qso_plot, marker='*', s=100, color='k')
    plotting.end_plot('plots/on_sky_distribution.png')
