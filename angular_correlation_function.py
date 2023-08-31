"""Angular two-point correlation function as defined in Ota, 2018, section 4."""

import numpy as np
import pylab as plt
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import angular_separation
from regionfy_catalog import load_region, get_region_mask
import plotting

def count_differences(ra_array, dec_array):
    """counts the angular separation between every point"""
    seps = [angular_separation(ra_array[i], dec_array[i], ra_array[i+1:], dec_array[i+1:]) for i in range(len(ra_array) -1)]
    return np.concatenate(seps)

def count_differnces_diff_arrays(ra_array_1, dec_array_1, ra_array_2, dec_array_2):
    """gets the angular separation between the two arrays."""
    seps = [angular_separation(ra_array_1[i], dec_array_1[i], ra_array_2, dec_array_2) for i in range(len(ra_array_1))]
    return np.concatenate(seps)

def sample_angular_distribution(distribution, theta_min, theta_max):
    """Finds the number of values between theta_min and theta_max in angular distribution."""
    cut = np.where((distribution <= theta_max) & (distribution > theta_min))[0]
    return len(cut)

def number_pairs(len_array):
    """Famous n(n-1)/2 forumala to normalize by the number of pairs"""
    return (len_array*(len_array +1 ))/2

def angular_function(data_data, data_random, random_random, theta_min, theta_max):
    """Angular correlation function for values between theta_min and theta_max."""
    d_d = sample_angular_distribution(data_data, theta_min, theta_max)
    d_r = sample_angular_distribution(data_random, theta_min, theta_max)
    r_r = sample_angular_distribution(random_random, theta_min, theta_max)

    d_d_norm = d_d/len(data_data)
    d_r_norm = d_r/len(data_random)
    r_r_norm = r_r/len(random_random)

    acf = (d_d_norm - (2 * d_r_norm) + r_r_norm)/r_r_norm
    u_acf = (1 + acf)/np.sqrt(d_d)  # see equation five of ota 2018
    return acf, u_acf

def full_angular_separation_function(data_data, data_random, random_random, theta_bins):
    """works out the angular sepration function over given theta bins."""
    y_values = []
    average_x_values = []
    y_err_values = []
    for i in range(len(theta_bins) -1):
        y_val, y_err = angular_function(data_data, data_random, random_random, theta_bins[i], theta_bins[i+1])
        x_val = (theta_bins[i] + theta_bins[i+1]) / 2
        y_values.append(y_val)
        average_x_values.append(x_val)
        y_err_values.append(y_err)
    return average_x_values, y_values, y_err_values


if __name__ == '__main__':
    INFILE = 'candidates_e.txt'
    DECAM_REGION_FILE = 'DECAM.reg'
    FITS_FILE  = '../correct_stacks/N964/n964.fits'

    #INFILE = 'candidates_cdfs_e.txt'
    #DECAM_REGION_FILE = '../CDFS_LAGER/DECAM_CDFS_FULL.reg'
    #FITS_FILE = 'CDFS_MASK.fits'

    hdul = fits.open(FITS_FILE)
    wcs = WCS(hdul[0].header)
    ra_candidates, dec_candidates = np.loadtxt(INFILE, unpack=True)
    decam_region = load_region(DECAM_REGION_FILE)

    ra_min, ra_max = 355.95*u.deg.to(u.rad), 358.3*u.deg.to(u.rad)
    dec_min, dec_max = -31.82*u.deg.to(u.rad), -29.84*u.deg.to(u.rad)

    #ra_min, ra_max = 53.8853418 * u.deg.to(u.rad), 51.4000866 * u.deg.to(u.rad)
    #dec_min, dec_max = -29 * u.deg.to(u.rad), -27*u.deg.to(u.rad)


    """
    Populate the random galaxies within the DECcam region.
    100 000, keeping with Ota, 2018.
    """

    NUMBER_RANDOM_GALS=1000
    ra_rad = np.random.uniform(ra_min,ra_max, NUMBER_RANDOM_GALS)
    dec_rad = np.arcsin(np.random.uniform(np.sin(dec_min), np.sin(dec_max), NUMBER_RANDOM_GALS))
    ra_deg = ra_rad*u.rad.to(u.deg)
    dec_deg = dec_rad*u.rad.to(u.deg)
    ra_pix, dec_pix = wcs.world_to_pixel_values(ra_deg, dec_deg)

    msk = get_region_mask(ra_pix, dec_pix, decam_region)
    ra_random, dec_random = ra_deg[msk], dec_deg[msk]

    print('DD')
    data_data_dist = count_differences(ra_candidates, dec_candidates)
    print('DR')
    data_random_dist = count_differnces_diff_arrays(ra_candidates, dec_candidates, ra_random, dec_random)
    print('RR')
    random_random_dist = count_differences(ra_random, dec_random)


    x, y, yerr = full_angular_separation_function(data_data_dist, data_random_dist, random_random_dist, np.arange(0,1.9,0.1))


    plotting.start_plot(r'$\theta$ [deg]', r'$\omega$ ($\theta$)')
    plt.errorbar(x, y, yerr=yerr, fmt='o')
    plt.scatter(x, y)
    plt.axhline(0, ls ='--', color='k', alpha=0.5)
    plotting.end_plot('plots/angular_tpcf.png')
