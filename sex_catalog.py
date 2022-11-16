"""Main class for sextractor catalogs."""

import numpy as np
import pandas as pd
from find_stars import cross_match
import sex_utils as su

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
