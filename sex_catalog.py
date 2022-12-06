"""Main class for sextractor catalogs."""

import numpy as np
import pylab as plt
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
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

    def to_region_file(self, decam_file: str, outfile: str, color: str = 'orange', only_galaxies = False) -> None:
        """Converts the points into a region file."""
        hdu = fits.open(decam_file)
        wcs, _ = find_optimal_celestial_wcs(hdu[1:])

        ra_positions = self.catalog['ALPHAPEAK_J2000'].values
        dec_positions = self.catalog['DELTAPEAK_J2000'].values
        x_positions, y_positions = wcs.world_to_pixel_values(ra_positions, dec_positions)

        if only_galaxies is True:
            cut = np.where(self.catalog['CLASS_STAR'].values < 0.5)
            x_positions = x_positions[cut]
            y_positions = y_positions[cut]
        write_positions_as_region_file(np.transpose((x_positions, y_positions)), outfile, color)

    def remove_sources_based_on_exposure_map(self, exp_map: str) -> None:
        """Uses an expsore map to exlude sources with no exposure times."""
        hdul = fits.open(exp_map)
        exptime = np.zeros(len(self.catalog['ALPHAPEAK_J2000']))
        for hdu in hdul[1:]:
            wcs = WCS(hdu.header)
            x_pixels, y_pixels = wcs.world_to_pixel_values(
                self.catalog['ALPHAPEAK_J2000'], self.catalog['DELTAPEAK_J2000'])
            x_pixels, y_pixels = x_pixels.astype(int), y_pixels.astype(int)
            cut_y = np.where((y_pixels < hdu.shape[0]) & (y_pixels >= 0))[0]
            cut_x = np.where((x_pixels < hdu.shape[1]) & (x_pixels >= 0))[0]
            final_cut = np.intersect1d(cut_x, cut_y).astype(int)
            vals = [hdu.data[y_pixels[index], x_pixels[index]] for index in final_cut]
            exptime[final_cut] = vals
        self.catalog['exptime'] = exptime
        self.catalog = self.catalog[self.catalog.exptime != 0]


if __name__ == '__main__':
    DECAM_FITS_FILE = \
        '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/c4d_211021_003940_osj_i_vik1.fits.fz'
    SEX_FILE = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/test.cat'
    i_cat = SExtractorCat(SEX_FILE)
    i_cat.to_region_file(DECAM_FITS_FILE, 'delete_region_file.reg', 'red')
