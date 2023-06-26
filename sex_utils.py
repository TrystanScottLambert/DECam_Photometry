"""Module to handle sextracor inputs and outputs."""

from typing import List, Tuple
from pandas import DataFrame

def _get_column_names(read_line_object: List) -> List:
    """Reads the header info of the sextractor catalog."""
    header = [line.split()[2] for line in read_line_object if line[0] == '#']
    return header

def _get_rows(read_line_object: List) -> List:
    """Takes the readline object and reads in the rows."""
    data = [list(map(float, line.split())) for line in read_line_object if line[0] != '#']
    return data

def split_names_and_data(read_line_object: List) -> Tuple[List, List]:
    """Takes the read in sextractor file and splits headers and columns."""
    header = _get_column_names(read_line_object)
    data = _get_rows(read_line_object)
    return header, data

def read_cat(sextractor_catalog: str) -> DataFrame:
    """Reads in the sextractor catalog."""
    with open(sextractor_catalog, encoding='utf8') as file:
        lines = file.readlines()
    column_names, data = split_names_and_data(lines)
    data_frame = DataFrame(data, columns = column_names)
    return data_frame


if __name__ == '__main__':
    INFILE = '/home/trystan/Desktop/Work/PhD/DECAM/correct_stacks/i/test.cat'
    cat_df = read_cat(INFILE)
