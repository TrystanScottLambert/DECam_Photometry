"""
Script to create a masked fits file and area determination. 
Used for plotting the actual used area in the source detection.

Also used for determining the actual area that is being used for 
source detection.
"""

import numpy as np
from rich.progress import track
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from regions import PolygonPixelRegion, PixCoord

def get_weights_mask(file_name:str, weight: float) -> tuple[np.ndarray, np.ndarray]:
    """Returns where there cannot be any detection."""
    hdu = fits.open(file_name)
    data = hdu[0].data
    return data < weight

def clean_weights_file(file_name: str, weight: float) -> np.ndarray:
    """Reads in and boolenfies the array"""
    hdu = fits.open(file_name)
    data = hdu[0].data
    data[data<weight] = 0
    data[data>weight] = 1
    return data

def calculate_vertices(center_x: float, center_y: float, width: float, height: float) -> tuple[list[float], list[float]]:
    """
    Calculate the vertices of a rectangle given its center coordinates, width, and height.

    Parameters:
        center_x (float): The x-coordinate of the center of the rectangle.
        center_y (float): The y-coordinate of the center of the rectangle.
        width (float): The width of the rectangle.
        height (float): The height of the rectangle.

    Returns:
        Tuple[List[float], List[float]]: A tuple containing two lists - x_vertices and y_vertices.
            x_vertices (List[float]): A list of x-coordinates of the four vertices of the rectangle.
            y_vertices (List[float]): A list of y-coordinates of the four vertices of the rectangle.
    """
    half_width = width / 2
    half_height = height / 2

    # Calculate the x and y coordinates of the four vertices
    x_vertices = [center_x - half_width, center_x + half_width, center_x + half_width, center_x - half_width]
    y_vertices = [center_y - half_height, center_y - half_height, center_y + half_height, center_y + half_height]

    return x_vertices, y_vertices

def parse_polygon_string(line: str) -> tuple[list, list]:
    """Reads in a region file string and returns the verticies as x, y for polygons"""
    coords = list(map(float, line.strip().split('(')[1].split(')')[0].split(',')))
    coords_2d = np.array(coords).reshape(-1, 2)
    x_coords = coords_2d[:, 0].astype(int)
    y_coords = coords_2d[:, 1].astype(int)
    return x_coords, y_coords

def parse_box_string(line: str) -> tuple[list, list]:
    """Reads in a region file string and returns the verticies as x, y for boxes"""
    values = np.array(line.split('(')[-1].split(')')[0].split(',')).astype(float)
    x_coords, y_coords = calculate_vertices(values[0], values[1], values[2], values[3])
    return x_coords, y_coords

def read_region_file(file_name: str) -> list[PolygonPixelRegion]:
    """Reads the region file and returns PolygonPixelRegion objects"""
    regions = []
    with open(file_name, encoding='utf8') as file:
        lines = file.readlines()

    for line in lines:
        if 'box' in line:
            x_coords, y_coords = parse_box_string(line)
            verticies = PixCoord(x = x_coords, y = y_coords)
            region = PolygonPixelRegion(vertices=verticies)
            regions.append(region)
        if 'polygon' in line:
            x_coords, y_coords = parse_polygon_string(line)
            verticies = PixCoord(x = x_coords, y = y_coords)
            region = PolygonPixelRegion(vertices=verticies)
            regions.append(region)

    return regions

if __name__ == '__main__':
    INFILE = '../CDFS_LAGER/cdfs_masks.reg'
    n964_file = '../CDFS_LAGER/n964_weight.fits'
    i_file = '../CDFS_LAGER/i_weight.fits'
    z_file = '../CDFS_LAGER/z_weight.fits'
    hdu = fits.open(n964_file)
    n964 = get_weights_mask(n964_file, 0.25)
    i = get_weights_mask(i_file, 0.4)
    z = get_weights_mask(z_file, 0.17)

    thing = np.ones(hdu[0].data.shape)
    thing[n964] = 0
    thing[i] = 0
    thing[z] = 0

    wcs = WCS(fits.open(n964_file)[0].header)
    regions  = read_region_file(INFILE)
    decam_region = read_region_file('../CDFS_LAGER/DECAM_CDFS.reg')[0]
    masks = [region.to_mask() for region in track(regions, 'making masks')]
    real_masks = [mask.to_image(hdu[0].data.shape) for mask in track(masks, 'to data shape')]
    region_mask = np.zeros(hdu[0].data.shape)
    for mask in real_masks:
        region_mask += mask
    region_mask = np.abs(region_mask -1)
    print('Made it here')
    final_mask = region_mask * thing

    print('Trying decam mask')
    decam_mask = decam_region.to_mask()
    decam_real_mask = decam_mask.to_image(hdu[0].data.shape)
    print('WE DID IT')
    
    print('Adding everything together')
    final_mask = final_mask * decam_real_mask
    print('JUST NEED TO PLOT')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=wcs)
    ax.imshow(final_mask)
    plt.show()

    hdu[0].data = final_mask
    hdu.writeto('CDFS_MASK.fits', overwrite=True)