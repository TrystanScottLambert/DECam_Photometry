"""
This script reduces sex catalogs to only include sources within a specified region.
It uses the Pyregion package to define the region of interest based on a DS9 region file,
and applies the mask to the input catalogs using numpy arrays.
The edited catalogs are then overwritten in place.

Dependencies:
- Pyregion
- numpy

Usage:
- Ensure that the Pyregion and numpy packages are installed.
- Define the path to the DS9 region file and the input catalogs.
- Run the script.
"""

import glob
from regions import PixCoord, PolygonPixelRegion
import numpy as np


def load_region(file_name: str):
    """
    Reads in a ds9 region file and returns a PolygonPixelRegion object.

    Parameters
    ----------
    file_name : str
        The name of the ds9 region file to load.

    Returns
    -------
    region_pix : `~regions.PolygonPixelRegion`
        The PolygonPixelRegion object representing the loaded region.

    """
    vertices = parse_ds9_region_file(file_name)
    x_pix, y_pix = split_xy(vertices)
    vertices = PixCoord(x=x_pix, y=y_pix)
    region_pix = PolygonPixelRegion(vertices=vertices)
    return region_pix

def parse_ds9_region_file(file_name:str) -> np.ndarray:
    """
    Returns an array of vertices from a ds9 region file.

    Parameters
    ----------
    file_name : str
        The name of the ds9 region file to load.

    Returns
    -------
    vertices_float : numpy.ndarray
        The array of vertices as floats.

    """
    with open(file_name, encoding='utf8') as file:
        lines = file.readlines()
    verticies_string = lines[-1].split('(')[-1].split(')')[0].split(',')
    verticies_float = np.array(verticies_string).astype(float)
    return verticies_float

def split_xy(verticies: np.ndarray) -> tuple:
    """
    Splits an array of vertices into separate x and y arrays.

    Parameters
    ----------
    vertices : numpy.ndarray
        The array of vertices.

    Returns
    -------
    tuple
        A tuple of two numpy.ndarray objects, the x and y arrays.

    """
    x_pix,y_pix = [],[]
    for i, _ in enumerate(verticies):
        if i%2 == 0:
            x_pix.append(verticies[i])
        elif i%2 == 1:
            y_pix.append(verticies[i])
    return x_pix, y_pix


def is_pixel_within_region(x_coord, y_coord, region):
    """
    Returns True if the given pixel is inside the specified region.

    Args:
        x_coord (float): The x-coordinate of the pixel.
        y_coord (float): The y-coordinate of the pixel.
        region (PolygonPixelRegion): The region to test against.

    Returns:
        bool: True if the pixel is inside the region, False otherwise.
    """
    return PixCoord(x_coord, y_coord) in region

def get_region_mask(x_array, y_array, region):
    """
    Returns a boolean mask indicating which pixels are inside the specified region.

    Args:
        x_array (numpy.ndarray): The x-coordinates of the pixels.
        y_array (numpy.ndarray): The y-coordinates of the pixels.
        region (PolygonPixelRegion): The region to test against.

    Returns:
        numpy.ndarray: A boolean mask indicating which pixels are inside the region.
    """
    mask = [is_pixel_within_region(x_array[i], y_array[i], region) for i in range(len(x_array))]
    return mask

def get_header_body(input_file: str):
    """
    Separates the header and body of an ASCII file.

    Args:
        input_file (str): The path to the input file.

    Returns:
        Tuple[list[str], list[str]]: A tuple containing the header lines and the body lines of the file.
    """
    with open(input_file, encoding='utf8') as file:
        lines = file.readlines()

    header = [line for line in lines if line.startswith('#')]
    body = [line for line in lines if not line.startswith('#')]
    return header, body

class Catalog:
    """A class to read, edit, and write sex catalogs based on a specified region.

    Attributes:
    -----------
    input_file : str
        Path to the input catalog file.
    header : list
        List of header lines in the input catalog.
    body : list
        List of body lines in the input catalog.
    x_array : numpy.ndarray
        Array of x-coordinates from the input catalog.
    y_array : numpy.ndarray
        Array of y-coordinates from the input catalog.
    region : regions.PolygonPixelRegion
        Region of interest, defined by a DS9 region file and loaded as a PolygonPixelRegion object using Pyregion.
    mask : list of bool
        List of True/False values indicating whether each source in the catalog falls within the specified region.

    Methods:
    --------
    reduce_body():
        Removes rows from the body of the catalog that correspond to sources outside of the specified region.
    write_catalog():
        Overwrites the input catalog with the reduced version.
    """
    def __init__(self, input_file:str, region: str):
        """
        Parameters:
        -----------
        input_file : str
            Path to the input catalog file.
        region : str
            Path to the DS9 region file defining the region of interest.
        """
        self.input_file = input_file
        self.header, self.body = get_header_body(input_file)
        self.x_array, self.y_array = np.loadtxt(self.input_file, usecols=(2,3), unpack=True)
        self.region = load_region(region)
        self.mask = get_region_mask(self.x_array, self.y_array, self.region)
        self.reduce_body()
        self.write_catalog()

    def reduce_body(self):
        """
        Updates the body of the catalog to include only those rows corresponding
        to sources within the region of interest.

        Returns:
        None
        """
        self.body = list(np.array(self.body)[self.mask])

    def write_catalog(self):
        """
        Writes the edited catalog to the original input file, overwriting its
        contents. The file is opened in 'w' mode, so any previous contents of
        the file will be erased.

        Returns:
        None
        """
        with open(self.input_file, 'w', encoding='utf8') as file:
            for line in self.header:
                file.write(line)

            for line in self.body:
                file.write(line)


if __name__ == '__main__':
    REGION_FILE = 'DECAM.reg'
    infiles = glob.glob('../correct_stacks/N964/*.cat')

    for file in infiles:
        Catalog(file, REGION_FILE)
