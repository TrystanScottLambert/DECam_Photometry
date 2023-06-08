"""
Script which reads in the weights files and uses them to generate a mask. 
"""

from astropy.io import fits
import numpy as np


def read_in_weights(infile: str) -> np.ndarray:
    """Reads in the data from the weights files"""
    hdu = fits.open(infile)
    data = hdu[0].data
    return data

def mask_data(data: np.ndarray, limit) -> np.ndarray:
    """Creates a mask for the data."""
    msk = np.where(data < limit)
    data[msk] = 0
    return data

def create_weight_mask(weights: list[tuple]) -> np.ndarray:
    """Creates the master weight mask which can be used for reducing sexcats"""
    datum = [read_in_weights(weight[0]) for weight in weights]
    masked_datum = [mask_data(data, weights[i][1]) for i, data in enumerate(datum)]
    return masked_datum[0] * masked_datum[1] * masked_datum[2]

def mask_catalog(x_array: np.array, y_array: np.array, mask: np.ndarray) -> np.ndarray:
    """Determines which indicies in the given x and y arrays are masked."""
    cat_mask = [
        mask[int(y_array[i]) -1, int(x_val) -1] != 0 for i, x_val in enumerate(x_array)
        ]
    return cat_mask

def get_header_body(infile: str):
    """Splits the header and the body of text file."""
    with open(infile, encoding='utf8') as file:
        lines = file.readlines()

    header = [line for line in lines if line.startswith('#')]
    body = [line for line in lines if not line.startswith('#')]
    return header, body


class Catalog:
    """Stores and edits catalog"""
    def __init__(self, sex_catalog:str, weights: list[tuple]):
        self.sex_catalog = sex_catalog
        self.master_mask = create_weight_mask(weights)
        print('Number of good pixels is ', len(np.where(self.master_mask !=0)[0]))
        self.header, self.body = get_header_body(sex_catalog)
        self.x_array, self.y_array = np.loadtxt(self.sex_catalog, usecols=(2,3), unpack=True)
        self.mask = mask_catalog(self.x_array, self.y_array, self.master_mask)
        self.reduce_body()
        self.write_catalog()

    def reduce_body(self):
        """Only selects rows that aren't masked"""
        self.body = list(np.array(self.body)[self.mask])

    def write_catalog(self):
        """Writes the current version of the catalog to file."""
        with open(self.sex_catalog, 'w', encoding='utf8') as file:
            for line in self.header:
                file.write(line)

            for line in self.body:
                file.write(line)


if __name__ == '__main__':
    WEIGHTS = [
        ('../correct_stacks/N964/n964_weight.fits', 0.008),
        ('../correct_stacks/N964/i_weight.fits', 0.0008),
        ('../correct_stacks/N964/z_weight.fits', 0.001)]

    SEX_CATALOGS = ['../correct_stacks/N964/n964.cat',
               '../correct_stacks/N964/n964_135.cat',
               '../correct_stacks/N964/i.cat',
               '../correct_stacks/N964/z.cat']

    for sex_cat in SEX_CATALOGS:
        Catalog(sex_cat, WEIGHTS)