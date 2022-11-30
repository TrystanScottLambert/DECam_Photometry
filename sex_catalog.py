"""Main class for sextractor catalogs."""

import numpy as np
import pylab as plt
import pandas as pd
from astropy.io import fits
from reproject.mosaicking import find_optimal_celestial_wcs
from find_stars import cross_match
import sex_utils as su
from decam import write_positions_as_region_file

class SExtractorCat:
    """Main class for sextractor catalogs."""
    def __init__(self, sextractor_file_name: str):
        """Initializing class"""
        self.catalog = su.read_cat(sextractor_file_name)

    def cross_match_with_panstars(self, pan_stars_file_name:str):
        """reducing the sextractor catalog to only those which are in panstars."""
        panstars_catalog = pd.read_csv(pan_stars_file_name)
        ra_panstars = np.array(list(panstars_catalog['raMean']))
        dec_panstars = np.array(list(panstars_catalog['decMean']))

        ra_decam = np.array(list(self.catalog['ALPHAPEAK_J2000']))
        dec_decam = np.array(list(self.catalog['DELTAPEAK_J2000']))

        idx_decam, idx_panstars = cross_match(ra_decam, dec_decam, ra_panstars, dec_panstars)

        matched_panstars_catalog = panstars_catalog.iloc[idx_panstars]
        matched_decam_catalog = self.catalog[idx_decam]
        return matched_decam_catalog, matched_panstars_catalog

    def quick_look(self) -> None:
        """Plots the ra and dec scatter plot to make sure things look ok."""
        ra_decam = np.array(list(self.catalog['ALPHAPEAK_J2000']))
        dec_decam = np.array(list(self.catalog['DELTAPEAK_J2000']))
        plt.scatter(ra_decam, dec_decam)
        plt.show()

    def to_region_file(self, decam_file: str, outfile: str, color: str = 'orange') -> None:
        """Converts the points into a region file."""
        hdu = fits.open(decam_file)
        wcs, _ = find_optimal_celestial_wcs(hdu[1:])

        ra_positions = self.catalog['ALPHAPEAK_J2000'].values
        dec_positions = self.catalog['DELTAPEAK_J2000'].values
        x_positions, y_positions = wcs.world_to_pixel_values(ra_positions, dec_positions)
        write_positions_as_region_file(np.transpose((x_positions, y_positions)), outfile, color)

if __name__ == '__main__':
    DECAM_FITS_FILE = \
        '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/c4d_211021_003940_osj_i_vik1.fits.fz'
    SEX_FILE = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/test.cat'
    i_cat = SExtractorCat(SEX_FILE)
    i_cat.to_region_file(DECAM_FITS_FILE, 'delete_region_file.reg', 'red')
