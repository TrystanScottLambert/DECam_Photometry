"""
Module for fitting exponential/splines to SNR-type data.

This is essential for determining things like depth using 
Sextractor uncertainties.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from sex_catalog import SExtractorCat
from zero_points import zero_points


def exponential_func(x, a, b):
    """Exponential function for fitting."""
    return a * np.exp(-b * x)

def calculate_snr(mag_err: float) -> float:
    """Converts the magnitude error into a snr value."""
    return (2.5/np.log(10))/mag_err


INFILE = '../correct_stacks/N964/n964.cat'
n964_cat_all = SExtractorCat(INFILE).catalog
n964_cat = n964_cat_all[n964_cat_all['MAG_AUTO'] != 99]
mag, err = n964_cat['MAG_AUTO'].values, n964_cat['MAGERR_AUTO'].values

mag += zero_points.n964_band.mag_correct(1)
mag_sort = np.argsort(mag)
mag = mag[mag_sort]
err = err[mag_sort]
snr = calculate_snr(err)

invalid_mask = np.logical_or(np.isnan(err), np.isinf(err))
valid_indices = np.where(~invalid_mask)
mag = mag[valid_indices]
err = err[valid_indices]
snr = snr[valid_indices]


params, cov_matrix = curve_fit(exponential_func, xdata=mag, ydata=err)

# Extract the fitting parameters
a_fit, b_fit = params
a_err, b_err = np.sqrt(np.diag(cov_matrix))

# Generating finer x values for smoother plot
x_finer = np.linspace(mag[0], mag[-1], 1000)

# Evaluate the fitted exponential curve at the finer x values
y_fitted = exponential_func(x_finer, a_fit, b_fit)
y_upper = exponential_func(x_finer, a_fit + 5 *a_err, b_fit + 5 * b_err)
y_lower = exponential_func(x_finer, a_fit - 5* a_err, b_fit - 5 * b_err)


# Plotting the original data and the fitted exponential curve
plt.scatter(mag, err)
plt.plot(x_finer, y_fitted, label='Exponential Fit', color='red')
plt.fill_between(x_finer, y_lower, y_upper, color='red', alpha=0.2, label='Uncertainty')
plt.xlabel('Magnitude')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
plt.show()
