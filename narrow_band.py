"""Narrow band analysis and fitting."""

import pylab as plt
import pandas as pd
from sex_catalog import SExtractorCat


if __name__ == '__main__':
    INFILE_SEX = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/N964/test.cat'
    INFILE_Z = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_z.csv'
    INFILE_Y = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_y.csv'

    decam_all = SExtractorCat(INFILE_SEX)
    decam_catalog_z, pan_cat_z = decam_all.cross_match_with_panstars(INFILE_Z)
    decam_catalog_y, pan_cat_y = decam_all.cross_match_with_panstars(INFILE_Y)
    #decam_catalog_final = pd.merge(decam_catalog_y, decam_catalog_z)
    #panstars_catalog_final = pd.merge(pan_cat_y, pan_cat_z)
    y = decam_catalog_z['MAG_BEST'].values# - pan_cat_z['zMeanPSFMag'].values
    x = pan_cat_z['zMeanPSFMag'].values - pan_cat_z['yMeanPSFMag'].values
    plt.scatter(x, y)
    plt.show()
