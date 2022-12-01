"""Narrow band analysis and fitting."""

import pylab as plt
import numpy as np
from scipy import odr
from sex_catalog import SExtractorCat
import plotting


def straight_line(parameters, x_val):
    """straight line for fitting to data"""
    alpha, zpt = parameters
    return x_val * alpha + zpt

def fwhm_plot(sex_file: str):
    """makes the fwhm showing the difference fwhm."""
    sex_cat = SExtractorCat(sex_file)
    plotting.start_plot('FWHM', 'Counts')
    plt.hist(sex_cat.catalog['FWHM_IMAGE'], bins = 1000, color='r', lw=3, histtype='step')
    plt.hist(sex_cat.catalog['FWHMPSF_IMAGE'], bins = 1000, color='b', lw=3, histtype='step')
    plt.legend()
    print('FWHMPSF: ', np.median(sex_cat.catalog['FWHMPSF_IMAGE']))
    print('FWHM:', np.median(sex_cat.catalog['FWHM_IMAGE']))
    plotting.end_plot('fwhm.png')

def make_region_file(decam_file:str, sex_file: str) -> None:
    """Makes the region file quickly."""
    sex_cat = SExtractorCat(sex_file)
    sex_cat.to_region_file(decam_file,'N964.reg')

if __name__ == '__main__':
    INFILE_SEX = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/N964/test.cat'
    INFILE_Z = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_z.csv'
    INFILE_Y = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_y.csv'

    decam_all = SExtractorCat(INFILE_SEX)
    decam_catalog_z, pan_cat_z = decam_all.cross_match_with_panstars(INFILE_Z)
    decam_catalog_y, pan_cat_y = decam_all.cross_match_with_panstars(INFILE_Y)

    z = pan_cat_z['zMeanPSFMag'].values
    y = pan_cat_z['yMeanPSFMag'].values
    z_err = pan_cat_z['zMeanPSFMagErr'].values
    y_err = pan_cat_z['zMeanPSFMagErr'].values
    n964 = decam_catalog_z['MAG_BEST'].values
    n964_err = decam_catalog_z['MAGERR_BEST'].values
    cut = np.where(y != -999)[0]
    z, y, n964, z_err, y_err, n964_err = z[cut], y[cut], n964[cut],\
         z_err[cut], y_err[cut], n964_err[cut]

    delta_pan = z - y
    delta_decam = z - n964
    delta_pan_err = np.hypot(z_err, y_err)
    delta_decam_err = np.hypot(z_err, n964_err)
    cut_x = np.where((delta_pan > -0.2) & (delta_pan < 0.6))[0]
    cut_y = np.where((delta_decam > 28.679) & (delta_decam < 29.752))
    final_cut = np.intersect1d(cut_x, cut_y)

    straight_line_model = odr.Model(straight_line)
    data = odr.RealData(
        delta_pan[final_cut], delta_decam[final_cut],
        sx=delta_pan[final_cut], sy=delta_decam[final_cut])

    _odr = odr.ODR(data, straight_line_model, beta0=[0.81, 28.9])
    out = _odr.run()
    popt = out.beta
    perr = out.sd_beta
    print("fit parameter 1-sigma error")
    print("———————————–")
    for i, pop in enumerate(popt):
        print(f'{pop} +- {perr[i]}')

    NSTD = 5.
    popt_up = popt + NSTD * perr
    popt_dw = popt - NSTD * perr
    x_fit = np.linspace(-0.3, 0.7, 100)
    fit = straight_line(popt, x_fit)
    fit_up = straight_line(popt_up, x_fit)
    fit_dw= straight_line(popt_dw, x_fit)

    fig = plotting.start_plot(x_label = 'z - y', y_label = 'N964 - z')
    plt.errorbar(
        delta_pan[final_cut],
        delta_decam[final_cut],
        xerr=delta_pan_err[final_cut],
        yerr = delta_decam_err[final_cut],
        fmt = 'ko', alpha=0.3, capsize=0, zorder=0)

    plt.plot(x_fit, fit,'r', lw=2)
    plt.fill_between(x_fit, fit_up, fit_dw, alpha=.5, color='r')
    plt.xlim(-0.3, 0.7)
    plotting.end_plot('/home/trystan/Desktop/Work/PhD/main/plots/N964_zpt_calculation.png')

    N964_DIR = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/N964/'
    N964_SEX_CAT = N964_DIR + 'test.cat'
    N964_DECAM_CAT = N964_DIR + 'c4d_210831_050404_osj_N964_vik1.fits.fz'

    #fwhm_plot(N964_SEX_CAT)
    n964_sex_cat = SExtractorCat(N964_SEX_CAT)
    n964_sex_cat.to_region_file(N964_DECAM_CAT, N964_DIR + 'N964.reg')
