"""
Script to write the candidates for the cdfs and our decam data as csv files
so that they can be converted into latex tables.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy.units.quantity import Quantity
import astropy.units as u


REDSHIFT_QSO = 6.9
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)

ARCSEC_PER_KPC_COMOVING = COSMO.arcsec_per_kpc_comoving(REDSHIFT_QSO)
DEG_PER_MPC_COMOVING = ARCSEC_PER_KPC_COMOVING.to(u.deg / u.Mpc)

ARCSEC_PER_KPC_PROPER = COSMO.arcsec_per_kpc_proper(REDSHIFT_QSO)
DEG_PER_MPC_PROPER = ARCSEC_PER_KPC_PROPER.to(u.deg / u.Mpc)


RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24) * u.deg
DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1 *u.deg
qso_position = SkyCoord(ra = RA_QSO, dec = DEC_QSO)

def calcuate_distance_to_quasar(positions: SkyCoord) -> tuple[list[str], list[str], list[str]]:
    """Works out the distances to the quasar in angular, comoving, and physical"""
    angular_sep = qso_position.separation(positions) # degrees
    comoving_distance = angular_sep/DEG_PER_MPC_COMOVING
    proper_distance = angular_sep/DEG_PER_MPC_PROPER

    return _stringify(angular_sep), _stringify(comoving_distance), _stringify(proper_distance)

def _stringify(quantity_array: Quantity) -> list[str]:
    """Converts an array of quantities into a rounded list of strings"""
    string_list = [str(round(val.value,2)) for val in quantity_array]
    string_list_zeros = [string_val+'0'  if len(string_val) == 3 else string_val for string_val in string_list]
    return string_list_zeros

def _format_skycoord(skycoord_value: SkyCoord) -> tuple[str, str]:
    """Takes the SkyCoord object and returns the position in hh:mm:ss.s dd:mm:ss.s format"""
    ra = skycoord_value.ra.hms
    dec = skycoord_value.dec.dms
    ra_string = f'{int(ra.h)}:{int(ra.m)}:{round(ra.s,1)}'
    dec_string = f'{int(dec.d)}:{int(dec.m)}:{round(dec.s,1)}'
    return ra_string, dec_string


def load_positions(candidate_file_name: str) -> tuple[list[str], list[str], SkyCoord]:
    """Reads in the candidate file and returns a SkyCoord object and the positons in hms format"""
    ra, dec = np.loadtxt(candidate_file_name, unpack=True)
    c = SkyCoord(ra = ra*u.deg, dec = dec*u.deg)
    string_positions = [_format_skycoord(val) for val in c]
    ras = [string_pos[0] for string_pos in string_positions]
    decs = [string_pos[1] for string_pos in string_positions]
    return ras, decs, c

def make_file(infile: str, outfile: str) -> None:
    """Writes the file given value array containing"""
    ras, decs, c = load_positions(infile)
    ang_dist, co_dist, prop_dist = calcuate_distance_to_quasar(c)
    with open(outfile, 'w', encoding='utf8') as file:
        for ra, dec, ang, co, prop in zip(ras, decs, ang_dist, co_dist, prop_dist):
            file.write(f'{ra},{dec},{ang},{co},{prop} \n')


if __name__ == '__main__':
    US_FILE = 'candidates_e.txt'
    CDFS_FILE = 'candidates_cdfs_e.txt'

    make_file(US_FILE, 'candidates.csv')
    make_file(CDFS_FILE, 'candidates_cdfs.csv')
