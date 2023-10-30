"""Plotting the shapely LBG spectrum for the thesis."""


import numpy as np 
import pylab as plt
from synphot.reddening import etau_madau


from plotting import prettify_plot, end_plot


wavelength, flux = np.loadtxt('shapely_spectrum.txt', unpack=True)


fig = plt.figure(figsize = (2*3.54, 3.54), dpi = 600)
prettify_plot(r'Rest-Frame Wavelength [$\AA$]', r'$F_{\nu}$ [$\mu$Jy]')
plt.plot(wavelength, flux/1e-29)
end_plot('shapely_fig2.png')