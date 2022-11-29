"""Narrow band analysis and fitting."""

import pylab as plt
import numpy as np
from scipy.optimize import curve_fit
from sex_catalog import SExtractorCat


def straight_line(x_val, alpha, zpt):
    """straight line for fitting to data"""
    return x_val * alpha + zpt

if __name__ == '__main__':
    INFILE_SEX = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/N964/test.cat'
    INFILE_Z = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_z.csv'
    INFILE_Y = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_y.csv'

    decam_all = SExtractorCat(INFILE_SEX)
    decam_catalog_z, pan_cat_z = decam_all.cross_match_with_panstars(INFILE_Z)
    decam_catalog_y, pan_cat_y = decam_all.cross_match_with_panstars(INFILE_Y)
    #decam_catalog_final = pd.merge(decam_catalog_y, decam_catalog_z)
    #panstars_catalog_final = pd.merge(pan_cat_y, pan_cat_z)
    z = pan_cat_z['zMeanPSFMag'].values
    y = pan_cat_z['yMeanPSFMag'].values
    z_err = pan_cat_z['zMeanPSFMagErr'].values
    y_err = pan_cat_z['zMeanPSFMagErr'].values
    n964 = decam_catalog_z['MAG_BEST'].values
    n964_err = decam_catalog_z['MAGERR_BEST'].values
    cut = np.where(y != -999)[0]
    z, y, n964, z_err, y_err, n964_err = z[cut], y[cut], n964[cut], z_err[cut], y_err[cut], n964_err[cut]

    delta_pan = z - y
    delta_decam = z - n964
    delta_pan_err = np.hypot(z_err, y_err)
    delta_decam_err = np.hypot(z_err, n964_err)
    cut_x = np.where((delta_pan > -0.2) & (delta_pan < 0.6))[0]
    cut_y = np.where((delta_decam > 28.679) & (delta_decam < 29.752))
    final_cut = np.intersect1d(cut_x, cut_y)

    plt.errorbar(delta_pan, delta_decam, xerr=delta_pan_err, yerr = delta_decam_err, fmt = 'ko', alpha=0.3, capsize=0)
    plt.show()

    a_fit, cov = curve_fit(straight_line, delta_pan[final_cut], delta_decam[final_cut], sigma=delta_decam_err[final_cut], absolute_sigma=True)
    errs = np.sqrt(np.diag(cov))
    plt.errorbar(delta_pan[final_cut], delta_decam[final_cut], xerr=delta_pan_err[final_cut], yerr = delta_decam_err[final_cut], fmt = 'ko', alpha=0.3, capsize=0)
    xs = np.linspace(-0.2, 0.6, 100)
    ys = straight_line(xs, *a_fit)
    plt.plot(xs, ys, lw=3, color='r')
    plt.xlabel('z - y', fontsize=14)
    plt.ylabel('N964 - z', fontsize=14)
    plt.title(f'alpha: {round(a_fit[0],2)} +- {round(errs[1],4)}; zpt: {round(a_fit[1],2)} +- {round(errs[1],4)}')
    plt.show()
