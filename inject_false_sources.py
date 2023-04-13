"""
Module to determine the sample incompleteness of the decam images using,
similiar to the method used in section 4.1 of Hu et. al., (2019). 
Except we do not use the balrog package and rely purly on galsim.
Note that we are only inserting sources into the narrow band image and
the flux is assumed to be zero for the broad band images.
"""

from typing import List
import numpy as np
import galsim
from astropy.io import fits
from regions import PixCoord

from hu_2019_plot import FILTERS
from regionfy_catalog import load_region
from plot_onsky_distribution import REGION_FILE

SERSIC_INDEX = 1.5
HALF_LIGHT_RADIUS = 0.9 # From Hu et. al., (2019)
PSF = 1.2 # arcseconds
MIN_MAG = 18
MAX_MAG = 26


def get_filter_from_image_name(file_name:str) -> str:
    """Searches the name of the decam image to determine the filter."""
    return file_name.split('.fits')[0].split('/')[-1]

class DecamImage:
    """Main class for the decam images."""
    def __init__(self, file_name: str) -> None:
        self.infile = file_name
        self.hdul = fits.open(file_name)
        self.filter = get_filter_from_image_name(file_name)
        self.zpt = FILTERS[self.filter].zpt

    @property
    def pixel_scale(self) -> float:
        """Gets the pixel scale of the image from the header."""
        return np.abs(float(self.hdul[0].header['PC1_1']))*3600 # arcseconds / pixel

    def calculate_counts_from_mag(self, mag: float) -> float:
        """Determines the counts that would be in the image given a certain magnitude."""
        return 10**((mag - self.zpt)/(-2.5))

    def generate_lae(self, mag:float) -> np.ndarray:
        """creates a lyman-alpha emitter galaxy image which can be injected into an image."""
        flux = self.calculate_counts_from_mag(mag)
        gal = galsim.Sersic(SERSIC_INDEX, HALF_LIGHT_RADIUS, flux=flux)
        psf = galsim.Gaussian(flux=1., sigma=PSF)
        final = galsim.Convolve([gal, psf])
        return final.drawImage(scale=self.pixel_scale).array

    def generate_many_laes(self, number_of_laes: int) -> List[np.ndarray]:
        """
        Creates several different lyman-alpha emitters from a mag-dist from Hu et. al., (2019).
        """
        mags = np.random.uniform(MIN_MAG, MAX_MAG, number_of_laes) # Hu ranges from 21 to 26 in magnitude
        laes = [self.generate_lae(mag) for mag in mags]
        return laes, mags

    def insert_lae(self, lae_data: np.ndarray, xpix_pos: int, ypix_pos: int) -> None:
        """
        Takes the smaller 2d array and inserts it into the data field.
        x_pos, y_pos need to be the corner of the image NOT THE CENTER.
        """
        x_len, y_len = lae_data.shape
        self.hdul[0].data[ypix_pos:ypix_pos+y_len, xpix_pos: xpix_pos+x_len] += lae_data

    def _generate_positions_in_region(self, number_of_positions: int) -> List:
        """Randomly selects positions within a give region."""
        region = load_region(REGION_FILE)
        y_len, x_len = self.hdul[0].data.shape
        counter = 0
        final_x_values = []
        final_y_values = []
        while counter < number_of_positions:
            y_positon = np.random.randint(0, y_len)
            x_position = np.random.randint(0, x_len)
            if PixCoord(x_position, y_positon) in region:
                counter += 1
                final_x_values.append(x_position)
                final_y_values.append(y_positon)
        return final_x_values, final_y_values

    def populate_image_with_mock_lae(self, number_of_mocks: int) -> List[int]:
        """
        Populates the image with mock lae emitters.
        Returns the x-positions, y-positions, and mags for analysis.
        """
        x_positions, y_positions = self._generate_positions_in_region(number_of_mocks)
        laes, mags = self.generate_many_laes(number_of_mocks)
        offset = laes[0].shape[0] / 2
        xs, ys = [], []
        for i, lae in enumerate(laes):
            self.insert_lae(lae, x_positions[i], y_positions[i])
            xs.append(x_positions[i] + offset)
            ys.append(y_positions[i] + offset)
        return xs, ys, mags

    def write(self):
        """Writes the injected field as a fits image ready for sextractor."""
        self.hdul.writeto(self.infile.split('.fits')[0] + '.injected.fits', overwrite=True)


if __name__ == '__main__':
    INFILE = '../correct_stacks/N964/n964.fits'
    n_band = DecamImage(INFILE)
    xs, ys, magnitudes = n_band.populate_image_with_mock_lae(3000)
    n_band.write()
    with open('mock_lae_sources.txt', 'w', encoding='utf8') as file:
        file.write('#x_pix y_pix magnitude \n')
        for i, mag in enumerate(magnitudes):
            file.write(f'{xs[i]} {ys[i]} {mag} \n')
