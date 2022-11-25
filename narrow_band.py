"""Narrow band analysis and fitting."""

import pylab as plt
from sex_catalog import SExtractorCat

def n_964(filters, alpha, zpt):
    """convert n_964"""
    z_filter, y_filter = filters
    return z_filter - alpha*(z_filter - y_filter) - zpt

if __name__ == '__main__':
    INFILE_SEX = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/N964/test.cat'
    INFILE_Z = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_z.csv'
    INFILE_Y = '/home/trystan/Desktop/Work/PhD/PANSTARS/PANSTARS_y.csv'

    decam_all = SExtractorCat(INFILE_SEX)
    decam_catalog_z, pan_cat_z = decam_all.cross_match_with_panstars(INFILE_Z)
    decam_catalog_y, pan_cat_y = decam_all.cross_match_with_panstars(INFILE_Y)
    #decam_catalog_final = pd.merge(decam_catalog_y, decam_catalog_z)
    #panstars_catalog_final = pd.merge(pan_cat_y, pan_cat_z)
    delta_pan = pan_cat_z['zMeanPSFMag'].values - pan_cat_z['yMeanPSFMag'].values
    delta_decam = decam_catalog_z['MAG_BEST'].values - pan_cat_z['zMeanPSFMag'].values
    plt.scatter(delta_pan, delta_decam)
    plt.show()
