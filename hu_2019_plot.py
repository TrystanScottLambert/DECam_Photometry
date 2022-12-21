"""Python script to display the DECAM filters."""

import pylab as plt
import numpy as np
import matplotlib as mpl
from astropy.coordinates import SkyCoord
import astropy.units as u
import plotting

mpl.rcParams.update({'font.size': 2})
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10)
mpl.rc('ytick', labelsize=10)

def read_in_test_cats(file_name: str):
    """read in the test cats for DECCAM"""
    ra, dec, mag = np.loadtxt(file_name, unpack=True, usecols=(0, 1, 4))
    return ra, dec, mag

def remove_dumb_values(array_1, array_2):
    """Will remove any 99 values that exist"""
    cut_1 = np.where(array_1 != 99)[0]
    cut_2 = np.where(array_2 != 99)[0]
    cut = np.union1d(cut_1, cut_2)
    return cut

if __name__ == '__main__':
    FWHM = 9.5
    sigma = FWHM / 2.355
    MU = 964
    FUDGE = 1000

    wavelength, i_band, z_band = np.loadtxt('DECAM_filters.txt', unpack=True, usecols=(0, 4, 5))
    '''fig = plt.figure(figsize = (2 * 3.54, 3.54), dpi = 600)
    plt.plot(wavelength, i_band, lw=2, label = 'i')
    plt.plot(wavelength, z_band, lw=2, label = 'z')
    x_plot = np.linspace(600, wavelength[-1], 100000)
    plt.plot(x_plot, FUDGE * stats.norm.pdf(x_plot, MU, sigma), label='NB964', lw=2.1)
    #plt.axvline(959.29,ls = '--', color='k', lw = 3)
    plt.axvline(960.4,ls = '--', color='k')
    plt.xlim(600)
    plt.legend()
    plt.xlabel('Wavelength [nm]', fontsize = 12)
    plt.ylabel('Arbitary Units', fontsize = 12)
    plt.tick_params(which='both', width=2,direction='in')
    plt.tick_params(which='major', length=4, direction='in')
    plt.savefig('Filters.png', bbox_inches = 'tight', transparent = False)'''

    ##############################

    Z_FILE = '/media/trystan/TOSHIBA EXT/DECAM/correct_stacks/N964/z.cat'
    N_FILE = '/media/trystan/TOSHIBA EXT/DECAM/correct_stacks/N964/n964.cat'
    RA_qso = (23 + (48/60) + (33.34/3600)) * (360/24)
    DEC_qso = (30 + (54/60) + (10/3600)) * -1
    z_ra, z_dec, z_mag = read_in_test_cats(Z_FILE)
    n_ra, n_dec, n_mag = read_in_test_cats(N_FILE)
    z_mag += 30.538
    n_mag += 29.012
    #good_idx = remove_dumb_values(z_mag, n_mag)
     
    catalog = SkyCoord(ra = n_ra*u.deg, dec = n_dec*u.deg)
    c = SkyCoord(ra = RA_qso * u.deg, dec = DEC_qso * u.deg)
    idx, d2d, _ = c.match_to_catalog_sky(catalog)

    plotting.start_plot('N964 [Mag]', 'z - N964 [Mag]')
    plt.scatter(n_mag, z_mag - n_mag, s=1, color='k', alpha=0.5)
    plt.scatter(n_mag[idx], z_mag[idx] - n_mag[idx], marker='*', s=50, color='m')
    plt.xlim(16, 28) 
    plt.ylim(-1,12)
    plt.axhline(1.9, color='r', lw=1)
    plotting.end_plot('plots/hu_plot_z.png')