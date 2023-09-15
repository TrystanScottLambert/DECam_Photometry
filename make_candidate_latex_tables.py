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
    return [f'{val.value:.2f}' for val in quantity_array]


def _format_skycoord(skycoord_value: SkyCoord) -> tuple[str, str]:
    """Takes the SkyCoord object and returns the position in hh:mm:ss.s dd:mm:ss.s format"""
    ra = skycoord_value.ra.hms
    dec = skycoord_value.dec.dms
    ra_second_string = f"{ra.s:.2f}"
    ra_string = f'{str(int(ra.h)).zfill(2)}:{str(int(ra.m)).zfill(2)}:{ra_second_string.zfill(5)}'
    dec_second_string = f"{abs(dec.s):.2f}"
    dec_string = f'{str(int(abs(dec.d))).zfill(2)}:{str(abs(int(dec.m))).zfill(2)}:{dec_second_string.zfill(5)}'
    return ra_string, dec_string


def load_positions(candidate_file_name: str) -> tuple[list[str], list[str], SkyCoord]:
    """Reads in the candidate file and returns a SkyCoord object and the positons in hms format"""
    ra, dec = np.loadtxt(candidate_file_name, unpack=True)
    c = SkyCoord(ra = ra*u.deg, dec = dec*u.deg)
    string_positions = [_format_skycoord(val) for val in c]
    ras = [string_pos[0] for string_pos in string_positions]
    decs = [string_pos[1] for string_pos in string_positions]
    return ras, decs, c


def read_luminosity_file(file_name: str) -> tuple[list[str], list[str], list[str]]:
    """
    Reading in the luminosity files for the candidates.
    MAKE SURE TO HAVE RUN lyman_luminosities.py before using this.
    """
    mag, mag_err, lya_lums, lya_lums_errs, sfrs, sfrs_errs = np.loadtxt(file_name, unpack=True)
    mag, mag_err = np.around(mag, 2), np.around(mag_err, 2)
    mag_err_str = [f"{err:.2f}" for err in mag_err]
    mag_str = [f"{m:.2f}"for m in mag]
    l, le, s, se = np.around(lya_lums, 1), np.around(lya_lums_errs,1), np.around(sfrs, 1), np.around(sfrs_errs,1)
    mag_strings = [f'{m} $\pm$ {err}' for m, err in zip(mag_str, mag_err_str)]
    lya_strings = [f'{lum} $\pm$ {err}' for lum, err in zip(l, le)]
    sfr_strings = [f'{sfr} $\pm$ {err}' for sfr, err in zip(s, se)]
    return lya_strings, sfr_strings, mag_strings


def make_file(infile: str, outfile: str, lum_file:str=None) -> None:
    """Writes the file given value array containing"""
    ras, decs, c = load_positions(infile)
    ang_dist, co_dist, prop_dist = calcuate_distance_to_quasar(c)
    data = [np.arange(1, len(ras) +1),ras, decs, ang_dist, co_dist, prop_dist]

    if lum_file is not None:
        lya_string, sfr_string, mag_strings = read_luminosity_file(lum_file)
        data += [mag_strings, lya_string, sfr_string]
    
    data = np.array(data).T
    lines = []
    for dat in data:
        line =''
        for val in dat:
            line+=f'{val},'
        lines.append(line)
    
    with open(outfile, 'w', encoding='utf8') as file:
        for line in lines:
            file.write(line[:-1] +' \n')


if __name__ == '__main__':
    US_FILE = 'candidates_e.txt'
    CDFS_FILE = 'candidates_cdfs_e.txt'
    LUM_FILE_US = 'candidate_luminosities.dat'
    LUM_FILE_CDFS = 'candidate_luminosities_cdfs.dat'

    make_file(US_FILE, 'candidates.csv', LUM_FILE_US)
    make_file(CDFS_FILE, 'candidates_cdfs.csv', LUM_FILE_CDFS)
