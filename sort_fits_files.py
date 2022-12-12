""""Sort directory of data into SCIENCE BIAS AND FLAT fields.
This is being used specifically to organize the IMACS data for THELI"""

import os
import glob
from typing import List
from astropy.io import fits


FOLDERS = {
    'Bias': 'BIAS',
    'Flat': 'FLAT',
    'Object': 'SCIENCE',
}


def make_directories(directory: str) -> None:
    """Makes the SCIENCE, BIAS, and FLAT directores."""
    for folder in FOLDERS.values():
        try:
            os.mkdir(f'{directory}/{folder}')
        except FileExistsError:
            print(f'Folder "{folder}" already exist.')

def get_fits_type(fits_file: str) -> str:
    """Reads the header and returns the fits file type."""
    hdu = fits.open(fits_file)
    object_type = hdu[0].header['EXPTYPE']
    return object_type

def move_to_directory(image_name: str, directory: str) -> None:
    """Takes an image and moves into a given directory."""
    os.replace(f'{image_name}', f'{directory}/{image_name}')

def sort_file(image_name: str) -> None:
    """Sort the file into the correct folder."""
    file_type = get_fits_type(image_name)
    move_to_directory(image_name, FOLDERS[file_type])

def get_list_of_all_fits_files(directory: str) -> List[str]:
    """Returns a list of all the fits files in the directory."""
    return glob.glob(f'{directory}/*.fits')

def sort_directory(directory: str) -> None:
    """Sorts an entire directory."""
    make_directories(directory)
    files = get_list_of_all_fits_files(directory)
    for file in files:
        sort_file(file)
