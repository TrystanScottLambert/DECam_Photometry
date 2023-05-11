"""Plotting the onsky distribution and the distance from the qso distribution"""

import numpy as np
from matplotlib import pyplot as plt
from regions import PixCoord, CirclePixelRegion
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import angular_separation
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
decam_hdu = fits.open(DECAM_SHAPE_FILE)
decam_wcs = WCS(decam_hdu[0].header)

RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1
REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
ARCSEC_PER_KPC = COSMO.arcsec_per_kpc_proper(REDSHIFT_QSO)
DEG_PER_MPC = ARCSEC_PER_KPC.to(u.deg / u.Mpc)
DEG_PER_PIX = np.abs(decam_hdu[0].header['PC2_2'])
REGION_FILE = 'decam.reg'
region_decam_fov = load_region(REGION_FILE)



if __name__ == '__main__':
    INFILE = 'candidates.txt'
    ra, dec = np.loadtxt(INFILE, unpack=True)

    ra_plot, dec_plot = decam_wcs.world_to_pixel_values(ra, dec)
    ra_qso_plot, dec_qso_plot = decam_wcs.world_to_pixel_values(RA_QSO, DEC_QSO)

    #On sky distribution plot.
    region_1 = set_region(ra_qso_plot, dec_qso_plot, 1)
    region_10 = set_region(ra_qso_plot, dec_qso_plot, 10)
    region_20 = set_region(ra_qso_plot, dec_qso_plot, 20)

    fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
    ax = fig.add_subplot(projection = decam_wcs)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.imshow(decam_hdu[0].data, alpha=0)
    region_1.plot(ax=ax, color='red', lw=2.0, ls=':')
    region_10.plot(ax=ax, color='red', lw=2.0, ls=':')
    region_20.plot(ax=ax, color='red', lw=2.0, ls=':')
    region_decam_fov.plot(ax = ax, color='k', lw=2.0)
    ax.scatter(ra_plot, dec_plot)
    ax.scatter(ra_qso_plot, dec_qso_plot, marker='*', s=100, color='k')
    plotting.end_plot('plots/on_sky_distribution.png')

    #on sky surface distribution plot.
    radii = np.arange(0,21,1) # Mpc
    average_radii_mpc = np.array([(radii[i] + radii[i+1])/2 for i in range(len(radii) -1)])
    average_radii_arcsec = average_radii_mpc * DEG_PER_MPC

    areas = np.array([np.pi * (radii[i+1]**2 -radii[i]**2) for i in range(len(radii) -1)])
    areas_deg = areas * (DEG_PER_MPC**2)

    angular_seperation_deg = angular_separation(RA_QSO, DEC_QSO, ra, dec)
    angular_seperation_pMpc = angular_seperation_deg / DEG_PER_MPC

    counts, _ = np.histogram(angular_seperation_pMpc.value, bins = radii)
    null_count_values = np.where(counts==0)[0]
    non_null_count_values = np.where(counts!=0)[0]
    counts[null_count_values] = 2 # poisson distribution means a zero count is less than 2 counts.
    y = counts/areas
    y_err = np.sqrt(counts)/areas

    fig = plt.figure(figsize = (3.54, 3.54/2), dpi = 600)
    ax = fig.add_subplot(111)
    #ax.errorbar(average_radii_mpc, y, yerr = y_err, fmt='ok', ecolor='r', ms = 4, capsize=2)
    ax.errorbar(average_radii_mpc[non_null_count_values], y[non_null_count_values], yerr = y_err[non_null_count_values], fmt='ok', ecolor='r', ms = 2, capsize=2)
    ax.errorbar(average_radii_mpc[null_count_values], y[null_count_values], yerr = y_err[null_count_values], fmt='ok', ecolor='r', ms = 2, capsize=2, uplims=True)
    ax.set_xlabel('Distance from QSO [pMpc]')
    ax.set_ylabel(r'Surface Density [pMpc$^{-2}$]')
    ax.minorticks_on()
    ax.tick_params(which='both', width=1.2,direction='in')
    ax.tick_params(which='major', length=3, direction='in')
    ax.axvline(14.04, ls=':', color='k', alpha=0.3)
    Y_LIM_MPC = 0.7
    PER_DEGREE_CONVERSION_FACTOR=1./(DEG_PER_MPC.value**2)
    #ax.set_ylim(top=Y_LIM_MPC, bottom=-0.05)
    ax.set_yscale('log')

    ax1 = ax.twiny()
    ax1.errorbar(average_radii_arcsec, y, yerr=y_err, fmt='ok', alpha=0)
    ax1.set_xlabel('Distance from QSO [deg]')
    ax1.minorticks_on()
    ax1.tick_params(which='both', width=1.2,direction='in')
    ax1.tick_params(which='major', length=3, direction='in')

    ax2 = ax.twinx()
    ax2.errorbar(average_radii_mpc, counts/areas_deg.value, yerr=np.sqrt(counts)/areas_deg.value, alpha=0)
    ax2.set_ylabel(r'Surface Density [deg$^{-2}$]')
    ax2.minorticks_on()
    ax2.tick_params(which='both', width=1.2,direction='in')
    ax2.tick_params(which='major', length=3, direction='in')
    #ax2.set_ylim(-0.05*PER_DEGREE_CONVERSION_FACTOR, Y_LIM_MPC*PER_DEGREE_CONVERSION_FACTOR)
    ax2.set_yscale('log')
    #plt.show()
    plotting.end_plot('plots/surface_density_log.png')
