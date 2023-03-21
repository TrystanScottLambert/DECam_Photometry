"""Python script to display the DECAM filters."""

import pylab as plt
import numpy as np
import matplotlib as mpl
from astropy.coordinates import SkyCoord
import astropy.units as u
import plotting
import postage_stamps as ps

mpl.rcParams.update({'font.size': 2})
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
mpl.rc('xtick', labelsize=10)
mpl.rc('ytick', labelsize=10)


def read_in_cats(file_name: str):
    """read in the test cats for DECCAM"""
    ra, dec, mag, mag_err = np.loadtxt(file_name, unpack=True, usecols=(0, 1, 4, 5))
    non_detections = np.where(mag == 99)[0]
    detections = np.where(mag != 99)[0]
    
    mag[detections] = mag[detections] + MAG_ZPTS[FILES[file_name]]
    mag_err[non_detections] = (2.5/np.log(10)) / 2 #factor of two is just for aesthitics
    
    return ra, dec , mag, mag_err

MAG_ZPTS = {
    'z': 30.538,
    'n964': 28.97,
    'i': 30.870,
}

LIMITING_MAGS = {
    '../correct_stacks/N964/z.cat': 25.5,
    '../correct_stacks/N964/i.cat': 25.5,
    '../correct_stacks/N964/n964.cat': 25.0
}

FILES = {
    '../correct_stacks/N964/z.cat': 'z',
    '../correct_stacks/N964/n964.cat': 'n964',
    '../correct_stacks/N964/i.cat': 'i'
}

if __name__ == '__main__':
    Z_FILE = '../correct_stacks/N964/z.cat'
    N_FILE = '../correct_stacks/N964/n964.cat'
    I_FILE = '../correct_stacks/N964/i.cat'
    RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
    DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1
    CANDIDATE_FILE = 'candidates.txt'

    candidates_ra, candidates_dec = np.loadtxt(CANDIDATE_FILE, usecols=(1,2), unpack=True)

    i_ra, i_dec, i_mag, i_mag_err = read_in_cats(I_FILE)
    z_ra, z_dec, z_mag, z_mag_err = read_in_cats(Z_FILE)
    n_ra, n_dec, n_mag, n_mag_err = read_in_cats(N_FILE)

    catalog = SkyCoord(ra = n_ra*u.deg, dec = n_dec*u.deg)
    c = SkyCoord(ra = RA_QSO*u.deg, dec = DEC_QSO*u.deg)
    idx, d2d, _ = c.match_to_catalog_sky(catalog)

    c_candidates = SkyCoord(ra = candidates_ra*u.deg, dec = candidates_dec*u.deg)
    idx_candidates, d2d_candidates, _ = c_candidates.match_to_catalog_sky(catalog)

    non_detections_z = np.where(z_mag[idx_candidates] == 99)[0]
    non_detections_i = np.where(i_mag[idx_candidates] == 99)[0]


    # Making the Hu plots
    plotting.start_plot('N964 [Mag]', 'z - N964 [Mag]')
    plt.scatter(n_mag, z_mag - n_mag, s=1, color='k', alpha=0.5)
    plt.scatter(n_mag[idx], z_mag[idx] - n_mag[idx], marker='*', s=50, color='m', label='QSO')
    plt.errorbar(n_mag[idx_candidates], z_mag[idx_candidates] - n_mag[idx_candidates],
                  color='r', xerr=n_mag_err[idx_candidates],
                yerr=np.hypot(n_mag_err[idx_candidates], z_mag_err[idx_candidates]), fmt='o', ms=2, elinewidth=1, label = 'Candidates')

    plt.errorbar(n_mag[idx_candidates][non_detections_z], LIMITING_MAGS[Z_FILE] - n_mag[idx_candidates][non_detections_z],
                  color='r', xerr=n_mag_err[idx_candidates][non_detections_z],
                yerr=np.hypot(n_mag_err[idx_candidates][non_detections_z], z_mag_err[idx_candidates][non_detections_z]),
                fmt='o', ms=2, elinewidth=1, lolims=True)
    plt.xlim(19, 25)
    plt.ylim(-1.1,6)
    plt.axhline(1.9, color='r', lw=1, ls='--')
    plt.legend(fontsize=8, loc=2)
    plotting.end_plot('plots/hu_plot_z.png')

    plotting.start_plot('N964 [Mag]', 'i - N964 [Mag]')
    plt.scatter(n_mag, i_mag - n_mag, s=1, color='k', alpha=0.5)
    plt.scatter(n_mag[idx], LIMITING_MAGS[I_FILE] - n_mag[idx], marker='*', s=50, color='m', label = 'QSO')
    plt.errorbar(n_mag[idx_candidates], i_mag[idx_candidates] - n_mag[idx_candidates],
                  color='r', xerr=n_mag_err[idx_candidates],
                yerr=np.hypot(n_mag_err[idx_candidates], i_mag_err[idx_candidates]), fmt='o', ms=2, elinewidth=1, label='Candidates')

    plt.errorbar(n_mag[idx_candidates][non_detections_i], LIMITING_MAGS[I_FILE] - n_mag[idx_candidates][non_detections_i],
                  color='r', xerr=n_mag_err[idx_candidates][non_detections_i],
                yerr=np.hypot(n_mag_err[idx_candidates][non_detections_i], i_mag_err[idx_candidates][non_detections_i]),
                fmt='o', ms=2, elinewidth=1, lolims=True)

    #plt.xlim(12.5, 28)
    plt.xlim(19, 25)
    plt.ylim(-1.1,10)
    plt.axhline(0.8, color='r', lw=1, ls='--')
    plt.legend(fontsize=8, loc=2)
    plotting.end_plot('plots/hu_plot_i.png')
    

    # Interactive plot that we can click and see what the galaxies look like.
    '''fig, ax = plt.subplots()
    coll = ax.scatter(n_mag, z_mag - n_mag, s=1, color='k', alpha=0.5, picker = True)
    ax.scatter(n_mag[idx], z_mag[idx] - n_mag[idx], marker='*', s=100, color='m')
    plt.scatter(n_mag[idx_candidates], z_mag[idx_candidates] - n_mag[idx_candidates], s = 20, color='r')
    ax.set_xlim(12.5, 28)
    ax.set_ylim(-2.4, 10)

    plt.axhline(1.9, color='r', lw=1)

    def on_pick(event):
        """what to do when clicking the button"""
        print('Position: ', n_ra[event.ind][0], n_dec[event.ind][0])
        ps.show_stamps(n_ra[event.ind][0], n_dec[event.ind][0])


    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()'''

    # plotting the magnitude values of the quasar itself.
    # comparison to the theoretical values from synphot
    '''print(i_ra[idx], i_dec[idx])
    print(d2d)
    print(RA_QSO, DEC_QSO)
    print('i: ', i_mag[idx])
    print('z: ', z_mag[idx])
    print('n964: ',n_mag[idx])
    print('--------------')
    print('i - n964: ', i_mag[idx] - n_mag[idx])
    print('i - z: ', i_mag[idx] - z_mag[idx])
    print('z - n964: ', z_mag[idx] - n_mag[idx])'''
