"""Making the postage stamp plot seen in many papers.
See: Figure 4, Banados, et. al., 2013"""

import numpy as np
import pylab as plt
from astropy.io import fits
from astropy.wcs import WCS

PAD = 15 # pixels

def cut_postage_stamp(r_a: float, dec: float, image):
    """Cutting out a postage stamp centered on r_a, and dec (in decimal degrees)"""
    wcs = WCS(image[0].header)
    x_pix, y_pix  = wcs.world_to_pixel_values(r_a, dec)
    data = image[0].data[
        int(y_pix)-PAD:int(y_pix)+PAD, int(x_pix)-PAD:int(x_pix)+PAD]
    return data

IMAGES = (
    '../correct_stacks/N964/i.fits',
    '../correct_stacks/N964/z.fits',
    '../correct_stacks/N964/n964.fits',
)

FITS_OBJECTS = [fits.open(image) for image in IMAGES]

if __name__ == '__main__':
    INFILE = 'candidates.txt'
    ra, dec = np.loadtxt(INFILE, usecols=(1,2), unpack=True)

    fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
    gs = fig.add_gridspec(len(ra), 3, hspace=0, wspace=0, width_ratios=[1, 1, 1],
                           top=0.80, bottom=0.20, left=0.4, right=0.5)
    axes = gs.subplots(sharex='col', sharey='row')
    for i, r_a in enumerate(ra):
        for j in range(3):
            axes[i, j].imshow(cut_postage_stamp(r_a, dec[i], FITS_OBJECTS[j]), cmap='gray_r')
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])

    axes[0,0].set_title('i', fontsize=3)
    axes[0,1].set_title('z', fontsize=3)
    axes[0,2].set_title('N964',fontsize=3)

    plt.savefig('postage_stamps.png', transparent = False)
