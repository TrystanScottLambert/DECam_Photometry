"""Plotting the onsky distribution"""

import numpy as np
from matplotlib import pyplot as plt, patches
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
import plotting
from regions import PixCoord, CirclePixelRegion


def get_radius(physical_radius_mpc: float) -> float:
    """Given a physical radius in Mpc, return the radius in arcseconds"""
    arcsec_per_kpc = COSMO.arcsec_per_kpc_proper(REDSHIFT_QSO)
    deg_per_mpc = arcsec_per_kpc.to(u.deg / u.Mpc)
    physical_radius_deg = physical_radius_mpc * u.Mpc * deg_per_mpc
    return physical_radius_deg.value

DECAM_SHAPE_FILE = '../correct_stacks/N964/n964_weight.fits'
decam_hdu = fits.open(DECAM_SHAPE_FILE)
decam_wcs = WCS(decam_hdu[0].header)

RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1
REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
DEG_PER_PIX = np.abs(decam_hdu[0].header['PC2_2'])

INFILE = 'candidates.txt'
ra, dec = np.loadtxt(INFILE, unpack=True)

ra_plot, dec_plot = decam_wcs.world_to_pixel_values(ra, dec)
ra_qso_plot, dec_qso_plot = decam_wcs.world_to_pixel_values(RA_QSO, DEC_QSO)
region_1 = CirclePixelRegion(center=PixCoord(x=ra_qso_plot, y=dec_qso_plot), radius=get_radius(1) / DEG_PER_PIX)
region_10 = CirclePixelRegion(center=PixCoord(x=ra_qso_plot, y=dec_qso_plot), radius=get_radius(10) / DEG_PER_PIX)
region_20 = CirclePixelRegion(center=PixCoord(x=ra_qso_plot, y=dec_qso_plot), radius=get_radius(20) / DEG_PER_PIX)

fig = plt.figure()
ax = fig.add_subplot(projection = decam_wcs)
ax.imshow(decam_hdu[0].data, alpha=0)
region_1.plot(ax=ax, color='red', lw=2.0, ls=':')
region_10.plot(ax=ax, color='red', lw=2.0, ls=':')
region_20.plot(ax=ax, color='red', lw=2.0, ls=':')
ax.scatter(ra_plot, dec_plot)
ax.scatter(ra_qso_plot, dec_qso_plot, marker='*', s=100, color='k')

plt.show()