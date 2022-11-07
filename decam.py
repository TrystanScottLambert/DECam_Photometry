"""Main class for DECAM images"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import vstack
from reproject.mosaicking import find_optimal_celestial_wcs
from tqdm import tqdm
from find_stars import find_stars_single, cross_match

def write_positions_as_region_file(positions: list, outfile: str, color: str):
    """Takes a list of positions and puts them in the ds9 region format."""
    with open(outfile, 'w', encoding="utf-8") as file:
        for position in positions:
            file.write(f'point({position[0]},{position[1]}) # point=circle 20  color={color}\n')

class DecamImage:
    """Loads in DECAM image and produces a star catalog."""
    def __init__(self, file_name: str):
        self.hdu_list = fits.open(file_name)
        self.main_wcs, _ = find_optimal_celestial_wcs(self.hdu_list[1:])
        self.create_star_catalog()

    def create_star_catalog(self):
        """Finding the stars across the entire stacked image."""
        star_catalogs = [find_stars_single(hdu) for hdu in tqdm(self.hdu_list[1:])]
        self.star_catalog = vstack(star_catalogs)

    def _convert_irsf_to_pixels(self, ra_array, dec_array):
        """Takes RA and Dec positions and puts them into pixel coordinates of the image."""
        x_pos, y_pos = self.main_wcs.world_to_pixel_values(ra_array, dec_array)
        positions = np.transpose((x_pos, y_pos))
        return positions

    def create_decam_cat_regions_file(self, outfile: str):
        """Create decam star-cat as a region file."""
        ras = np.array(list(self.star_catalog['RA']))
        decs = np.array(list(self.star_catalog['Dec']))
        self.create_region_file(ras, decs, outfile)
    
    def create_region_file(self, ra_array, dec_array, outfile: str, color='green'):
        """Generalized function to write region files given ra and dec arrays."""
        positions = self._convert_irsf_to_pixels(ra_array, dec_array)
        write_positions_as_region_file(positions, outfile, color=color)

    def create_cat_regions_file(self, outfile: str):
        """Takes a star catalog and converts it into a .reg file for ds9."""
        ras = np.array(list(self.star_catalog['RA']))
        decs = np.array(list(self.star_catalog['Dec']))
        x_pos, y_pos = self.main_wcs.world_to_pixel_values(ras, decs)
        positions = np.transpose((x_pos, y_pos))
        with open(outfile, 'w', encoding="utf-8") as file:
            for position in positions:
                file.write(f'point({position[0]},{position[1]}) # point=circle 20 \n')

    def cross_match_with_panstars(self, pan_stars_catalog_file_name: str):
        """panstars catalog must be in correct form from prep_panstars.py"""
        panstars_catalog = pd.read_csv(pan_stars_catalog_file_name)
        ra_panstars = np.array(list(panstars_catalog['raMean']))
        dec_panstars = np.array(list(panstars_catalog['decMean']))
        idx_decam, idx_panstars = cross_match(self.star_catalog['RA'], self.star_catalog['Dec'],
                                              ra_panstars, dec_panstars)
        matched_panstars_catalog = panstars_catalog.iloc[idx_panstars]
        matched_decam_catalog = self.star_catalog[idx_decam]

        self.create_region_file(matched_panstars_catalog['raMean'],
                                matched_panstars_catalog['decMean'],
                                'panstars_matches.reg',
                                color = 'red')

        self.create_region_file(matched_decam_catalog['RA'],
                                matched_decam_catalog['Dec'], 'decam_matches.reg')
        return matched_decam_catalog, matched_panstars_catalog


if __name__ == '__main__':
    INFILE = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/c4d_211021_003940_osj_i_vik1_skysubtracted.fits.fz'
    PAN_CAT = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_i.csv'
    test = DecamImage(INFILE)
    test.create_decam_cat_regions_file('test.reg')
    thing = test.cross_match_with_panstars(PAN_CAT)
