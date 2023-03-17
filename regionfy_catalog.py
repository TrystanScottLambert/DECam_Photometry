""" Reducing sex catlogs to only include soruces within a region. """

from regions import PixCoord, PolygonPixelRegion
import numpy as np

def load_region(file_name: str):
    """reads in the ds9 region file as a Region object"""
    vertices = get_verticies(file_name)
    x_pix, y_pix = split_xy(vertices)
    vertices = PixCoord(x=x_pix, y=y_pix)
    region_pix = PolygonPixelRegion(vertices=vertices)
    return region_pix

def get_verticies(file_name:str):
    """Returns the verticies of ds9 region as an array"""
    file = open(file_name, encoding='utf8')
    lines = file.readlines()
    verticies_string = lines[-1].split('(')[-1].split(')')[0].split(',')
    verticies_float = np.array(verticies_string).astype(float)
    return verticies_float

def split_xy(verticy_array):
    """ Split the array of verticies into x and y arrays."""
    x_pix,y_pix = [],[]
    for i in range(len(verticy_array)):
        if i%2 == 0:
            x_pix.append(verticy_array[i])
        elif i%2 == 1:
            y_pix.append(verticy_array[i])
    return x_pix, y_pix


def is_in_region(x_coord, y_coord, region):
    """Tells you in the particular pixel is in the region"""
    return PixCoord(x_coord, y_coord) in region

def mask_in_region(x_array, y_array, region):
    """Returns a True/False array for if the positoins are in or out of the region"""
    mask = [is_in_region(x_array[i], y_array[i], region) for i in range(len(x_array))]
    return mask

def get_headerbody(infile:str):
    """Separates the header and body of an ascii file"""
    file = open(infile, encoding='utf8')
    lines = file.readlines()
    file.close()
    header = [line for line in lines if line[0][0] == '#']
    body = [line for line in lines if line[0][0] != '#']
    return header, body

class Catalog:
    """Class to read, edit, and write the catalogs"""
    def __init__(self, infile:str, region: str):
        self.infile = infile
        self.header, self.body = get_headerbody(infile)
        self.x_array, self.y_array = np.loadtxt(self.infile, usecols=(2,3), unpack=True)
        self.region = load_region(region)
        self.mask = mask_in_region(self.x_array, self.y_array, self.region)
        self.reduce_body()
        self.write_catalog()

    def reduce_body(self):
        """Only returns the rows which are not masked."""
        self.body = list(np.array(self.body)[self.mask])

    def write_catalog(self):
        """OVERWRITE the catalog"""
        file = open(self.infile, 'w', encoding='utf8')
        for line in self.header:
            file.write(line)

        for line in self.body:
            file.write(line)


if __name__ == '__main__':
    REGION_FILE = 'DECAM.reg'
    INFILE_N964 = '../correct_stacks/N964/n964.cat'
    INFILE_Z = '../correct_stacks/N964/z.cat'
    INFILE_I = '../correct_stacks/N964/i.cat'

    Catalog(INFILE_N964, REGION_FILE)
    Catalog(INFILE_Z, REGION_FILE)
    Catalog(INFILE_I, REGION_FILE)

