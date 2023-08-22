"""
Masking the catalogs using the created mask from make_mask_fits_file.
"""

import glob
from astropy.io import fits
import numpy as np

def get_header_body(input_file: str) -> tuple[list[str], list[str]]:
    """
    Separates the header and body of an ASCII file.

    Args:
        input_file (str): The path to the input file.

    Returns:
        Tuple[list[str], list[str]]: A tuple containing the header lines
        and the body lines of the file.
    """
    with open(input_file, encoding='utf8') as file:
        lines = file.readlines()

    header = [line for line in lines if line.startswith('#')]
    body = [line for line in lines if not line.startswith('#')]
    return header, body


class Mask:
    """
    Masks made by the maske_mask_fits_file.py These masks are fits files with 
    0 representing areas that are masked and 1s areas that are not.
    """
    def __init__(self, infile: str) -> None:
        """Initilizing"""
        self.hdu = fits.open(infile)
        self.data = self.hdu[0].data

    def is_pix_good(self, x_val: int, y_val: int) -> bool:
        """Determines if a given point is masked."""
        good = True
        if self.data[y_val-1, x_val-1] == 0:
            good = False
        return good

def mask_xy_vals(x_array: np.ndarray, y_array: np.ndarray, mask: Mask) -> list[bool]:
    """Loops through all values and determines if they are masked or not."""
    mask = [mask.is_pix_good(int(x_val), int(y_val)) for x_val, y_val in zip(x_array, y_array)]
    return mask


class Catalog:
    """Catalogs"""
    def __init__(self, infile: str) -> None:
        """Initilizing"""
        self.infile = infile
        self.header, self.body = get_header_body(self.infile)
        self.x_array, self.y_array = np.loadtxt(infile, usecols=(2,3), unpack=True)

    def mask_values(self, mask: Mask) -> None:
        """Masks the current values."""
        val_mask = mask_xy_vals(self.x_array, self.y_array, mask)
        self.body = list(np.array(self.body)[val_mask])

    def write_catalog(self) -> None:
        """
        Writes the edited catalog to the original input file, overwriting its
        contents. The file is opened in 'w' mode, so any previous contents of
        the file will be erased.

        Returns:
        None
        """
        with open(self.infile, 'w', encoding='utf8') as file:
            for line in self.header:
                file.write(line)

            for line in self.body:
                file.write(line)

if __name__ == '__main__':
    MASK_NAME = 'DECAM_MASK.fits'
    catalogs = glob.glob('../correct_stacks/N964/*.cat')
    cdfs_mask = Mask(MASK_NAME)

    for catalog in catalogs:
        cat = Catalog(catalog)
        cat.mask_values(cdfs_mask)
        cat.write_catalog()
