"""
Module which determines the completeness value based on mock source injection. 
See Figure 4, section 4.1, and f_{comp} of equation 5 in Hu et. al., (2019).
"""

import scipy
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import pylab as plt
from scipy.optimize import curve_fit
from inject_false_sources import MIN_MAG, MAX_MAG
from plotting import start_plot, end_plot

from zero_points import zero_points
from identify_candidates import calculate_snr, read_all
from depth_values import OUR_DEPTH
from decam import write_positions_as_region_file


def get_matches(target_catalog: SkyCoord, main_catalog: SkyCoord):
    """
    Finds the sources in the target catalog which have appropriate matches
    in the main_catalog. This is determined through a tolerance of 1".

    Returns the indicies of the main_catalog and then the indicies of the 
    target catalog.
    """
    idx, d2d, _ = target_catalog.match_to_catalog_sky(main_catalog)
    cut = np.where(d2d < 1 *u.arcsec)[0]
    return cut, idx[cut]


def remove_original_sources(recovered_catalog: SkyCoord, original_catalog: SkyCoord):
    """
    Matches the sources that are identical in the original catalog and removes them.
    If the positions are exactly the same as those in the original catalog then they have 
    to be sources from the original catalog.
    """

    _, d2d,_ = recovered_catalog.match_to_catalog_sky(original_catalog)
    cut = np.where(d2d != 0)[0]
    return cut

if __name__ == '__main__':
    MOCK_SOURCES_2_SEXCAT = '../correct_stacks/N964/n964_false.cat'
    ORIGINAL_SEXCAT = '../correct_stacks/N964/n964.cat'
    MOCK_SOURCES_FILE = 'mock_lae_sources.txt'
    FITS_FILE = '../correct_stacks/N964/n964.injected.fits'
    I_MOCK_CATALOG = '../correct_stacks/N964/i_false.cat'
    I_CATALOG = '../correct_stacks/N964/i.cat'
    Z_MOCK_CATALOG = '../correct_stacks/N964/z_false.cat'
    Z_CATALOG = '../correct_stacks/N964/z.cat'

    hdu = fits.open(FITS_FILE)
    wcs = WCS(hdu[0].header)
    x, y, mag = np.loadtxt(MOCK_SOURCES_FILE, unpack=True)
    flux = 10**((mag - zero_points.n964_band.mag_correct(1))/(-2.5))
    ra, dec = wcs.pixel_to_world_values(x, y)

    ra_original, dec_original, flux_original = np.loadtxt(ORIGINAL_SEXCAT, usecols=(0, 1, 6), unpack=True)
    ra_mock, dec_mock, x_mock, y_mock, mag_mock, flux_mock = np.loadtxt(MOCK_SOURCES_2_SEXCAT, usecols=(0, 1, 2, 3, 4, 6), unpack=True)
    mag_mock += zero_points.n964_band.mag_correct(1)

    original_cat = SkyCoord(ra = ra_original * u.deg, dec = dec_original * u.deg)
    recovered_cat = SkyCoord(ra = ra_mock * u.deg, dec = dec_mock * u.deg)
    cat = SkyCoord(ra = ra * u.deg, dec = dec * u.deg)

    # Match the injected catalog with the known mock positions.
    idx_recovered, idx_mock = get_matches(recovered_cat, cat)

    # Remove any sources that are actually from the original catalog
    possible_recovered_idx = remove_original_sources(recovered_cat[idx_recovered], original_cat)
    idx_found = idx_recovered[possible_recovered_idx]

    #Reading in the broadband data.
    mag_z, err_z, snr_z = read_all(Z_MOCK_CATALOG)
    mag_i, err_i, snr_i = read_all(I_MOCK_CATALOG)
    mag_n964, err_n964, snr_n964 = read_all(MOCK_SOURCES_2_SEXCAT)

    mag_z += zero_points.z_band.mag_correct(1)
    mag_i += zero_points.i_band.mag_correct(1)
    mag_n964 += zero_points.n964_band.mag_correct(1)

    #Preparing the values (including setting limiting magnitudes)
    mag_z, err_z, snr_z = mag_z[idx_found], err_z[idx_found], snr_z[idx_found]
    mag_i, err_i, snr_i = mag_i[idx_found], err_i[idx_found], snr_i[idx_found]
    mag_n964, err_n964, snr_n964 = mag_n964[idx_found], err_n964[idx_found], snr_n964[idx_found]


    mag_z = mag_n964 + 2.576
    mag_i[np.where(snr_i < 1)] = OUR_DEPTH.i_band.sigma_1
    mag_z[np.where(snr_z < 1)] = OUR_DEPTH.z_band.sigma_3
    


    #Narrowband Color selection
    color = mag_z - mag_n964
    significance = 2.5 * np.hypot(err_z, err_n964)
    first_cut = np.where(color > 0.78)[0]
    final_cut = []
    for idx in first_cut:
        if np.abs(color[idx]) > significance[idx]:
            final_cut.append(idx)
    narrow_band_cut = np.array(final_cut)


    #Continuum color selection
    color = mag_i - mag_z
    continuum_cut = np.where(color > 1.)[0]

    #i-band cut
    i_band_cut = np.where(snr_i < 2)[0]

    #Selection
    selection = np.intersect1d(i_band_cut,np.intersect1d(continuum_cut, narrow_band_cut))



    #Plotting


    def func(x, a, b, z, f):
        return a * scipy.special.erf(f*(x - z)) + b
    

   
    BINWIDTH = 0.2
    bins = np.arange(MIN_MAG, MAX_MAG + BINWIDTH, BINWIDTH)
    x_avg = [(bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]

    y_mock, _ = np.histogram(mag, bins=bins)
    y_recovered, _ = np.histogram(mag[idx_mock][possible_recovered_idx], bins=bins)

    y_completeness = y_recovered/y_mock
    popt, pcov = curve_fit(func, x_avg, y_completeness, maxfev=8000, p0=[-0.5,0.5,24,1])
    perr = np.sqrt(np.diag(pcov))
    popt_up = popt + perr
    popt_down = popt - perr

    print('popt', *popt)
    start_plot('NB964 Magnitude','Completeness')
    plt.errorbar(x_avg, y_completeness, yerr=np.sqrt(y_recovered)/y_mock, fmt='o', elinewidth=1, ms=2, color='k', ecolor='k')
    x_plotting = np.linspace(x_avg[0], x_avg[-1], 10000)
    plt.plot(x_plotting, func(x_plotting, *popt), 'r')
    plt.fill_between(x_plotting, func(x_plotting, *popt_down), func(x_plotting, *popt_up), color='r', alpha=0.3)
    end_plot('plots/completeness.png')
    plt.show()
