"""
Program for identifying stars in DECAM images.
"""

import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.table import Column
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from photutils.detection import DAOStarFinder

def cross_match(ra_new: np.array, dec_new: np.array, ra_catalog: np.array, dec_catalog: np.array):
    "Cross match to sets of coordinates."
    c_new = SkyCoord(ra = ra_new*u.deg, dec = dec_new*u.deg)
    c_catalog = SkyCoord(ra = ra_catalog*u.deg, dec = dec_catalog*u.deg)
    idx, d2d, _ = c_new.match_to_catalog_sky(c_catalog)
    max_separation = 1.0 * u.arcsec
    separation_constraint = d2d < max_separation
    c_idx = separation_constraint
    catalog_idx = idx[c_idx]
    return c_idx, catalog_idx


def find_stars_single(fits_hdu):
    "Finds stars in a fits data image and returns a catalog."
    _, median, std = sigma_clipped_stats(fits_hdu.data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=5, threshold=20*std)
    sources = daofind(fits_hdu.data - median)

    #have to get RA and Dec because x, y are local positions
    x_pos = np.array(list(sources['xcentroid']))
    y_pos = list(sources['ycentroid'])
    local_wcs = WCS(fits_hdu.header)
    ras, decs = local_wcs.pixel_to_world_values(x_pos, y_pos)
    col_ra = Column(name='RA', data = ras)
    col_dec = Column(name='Dec', data = decs)
    sources.add_columns([col_ra, col_dec])

    return sources
