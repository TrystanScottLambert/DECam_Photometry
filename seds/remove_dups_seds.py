"""
Script to remove any duplicate entries in any sed files

Specifically for plot_color_color_model.py
"""

import glob
import numpy as np

def remove_duplicates(file_name: str) -> None:
    """Opens the file, removes the duplicates, updates the files."""
    data = np.genfromtxt(file_name)
    _, unique_indices = np.unique(data[:, 0], return_index=True)
    unique_data = data[unique_indices]
    np.savetxt(file_name, unique_data, fmt='%f %f')

if __name__ == '__main__':
    files = glob.glob('*.sed')
    for file in files:
        remove_duplicates(file)
