"""Plotting the onsky distribution and the distance from the qso distribution"""

import numpy as np
from matplotlib import pyplot as plt
from regions import PixCoord, CirclePixelRegion
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
from regionfy_catalog import load_region
from astropy.visualization import ZScaleInterval
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


DECAM_SHAPE_FILE = '../CDFS_LAGER/CDFS_NB.fits.weight.fits'
decam_hdu = fits.open(DECAM_SHAPE_FILE)
decam_wcs = WCS(decam_hdu[0].header)

center_y, center_x = decam_hdu[0].shape
center_y, center_x = center_y/2, center_x/2
RA_QSO, DEC_QSO = decam_wcs.pixel_to_world_values(center_x, center_y)
REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
ARCSEC_PER_KPC = COSMO.arcsec_per_kpc_comoving(REDSHIFT_QSO)
DEG_PER_MPC = ARCSEC_PER_KPC.to(u.deg / u.Mpc)
DEG_PER_PIX = np.abs(decam_hdu[0].header['CD2_2'])
REGION_FILE = '../CDFS_LAGER/DECAM_CDFS_FULL.reg'
region_decam_fov = load_region(REGION_FILE)


if __name__ == '__main__':

    CDFS_MASK = fits.open('CDFS_MASK.fits')
    cdfs_mask_data = CDFS_MASK[0].data
    cdfs_mask_data[cdfs_mask_data == 0] = np.nan
    image = fits.open('../CDFS_LAGER/CDFS_NB.fits')
    masked_image = image[0].data * cdfs_mask_data
    zscale = ZScaleInterval()
    lower_lim, upper_lim = zscale.get_limits(masked_image)
    INFILE = 'candidates_cdfs_e.txt'
    ra, dec = np.loadtxt(INFILE, unpack=True)
    ra_snr, dec_snr = np.loadtxt('candidates_true_cdfs.txt', unpack=True)

    ra_plot, dec_plot = decam_wcs.world_to_pixel_values(ra, dec)
    ra_snr_plot, dec_snr_plot = decam_wcs.world_to_pixel_values(ra_snr, dec_snr)
    ra_qso_plot, dec_qso_plot = decam_wcs.world_to_pixel_values(RA_QSO, DEC_QSO)

    #On sky distribution plot.
    region_1 = set_region(ra_qso_plot, dec_qso_plot, 75)
    region_10 = set_region(ra_qso_plot, dec_qso_plot, 10)
    region_20 = set_region(ra_qso_plot, dec_qso_plot, 20)

    fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
    ax = fig.add_subplot(projection = decam_wcs)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    #ax.imshow(cdfs_mask_data, alpha=0.4, cmap='gray_r')
    ax.imshow(masked_image, cmap='gray_r', vmin=lower_lim, vmax=upper_lim)
    region_1.plot(ax=ax, color='red', lw=2.0, ls=':')
    #region_10.plot(ax=ax, color='red', lw=2.0, ls=':')
    #region_20.plot(ax=ax, color='red', lw=2.0, ls=':')
    region_decam_fov.plot(ax = ax, color='k', lw=2.0)
    ax.scatter(ra_plot, dec_plot, s=10, label = r'i$_{snr}$ > 2')
    #ax.scatter(ra_snr_plot, dec_snr_plot, s=10, color='r', marker='s', label=r'i$_{snr}$ < 2')
    ax.legend(frameon=False, loc=1)
    plotting.end_plot('plots/on_sky_distribution_cdfs.png')
