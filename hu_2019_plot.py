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




def read_in_test_cats(file_name: str):
    """read in the test cats for DECCAM"""
    ra, dec, mag = np.loadtxt(file_name, unpack=True, usecols=(0, 1, 4))
    return ra, dec, mag


MAG_ZPTS = {
    'z': 30.538,
    'n964': 29.012,
    'i': 30.870,
}

if __name__ == '__main__':
    Z_FILE = '../correct_stacks/N964/z.cat'
    N_FILE = '../correct_stacks/N964/n964.cat'
    I_FILE = '../correct_stacks/N964/i.cat'
    RA_QSO = (23 + (48/60) + (33.34/3600)) * (360/24)
    DEC_QSO = (30 + (54/60) + (10.0/3600)) * -1

    i_ra, i_dec, i_mag = read_in_test_cats(I_FILE)
    z_ra, z_dec, z_mag = read_in_test_cats(Z_FILE)
    n_ra, n_dec, n_mag = read_in_test_cats(N_FILE)

    z_mag += MAG_ZPTS['z']
    n_mag += MAG_ZPTS['n964']
    i_mag += MAG_ZPTS['i']

    catalog = SkyCoord(ra = n_ra*u.deg, dec = n_dec*u.deg)
    c = SkyCoord(ra = RA_QSO * u.deg, dec = DEC_QSO * u.deg)
    idx, d2d, _ = c.match_to_catalog_sky(catalog)

    print(i_ra[idx], i_dec[idx])
    print(d2d)
    print(RA_QSO, DEC_QSO)
    print('i: ', i_mag[idx])
    print('z: ', z_mag[idx])
    print('n964: ',n_mag[idx])
    print('--------------')
    print('i - n964: ', i_mag[idx] - n_mag[idx])
    print('i - z: ', i_mag[idx] - z_mag[idx])
    print('z - n964: ', z_mag[idx] - n_mag[idx])

        
    # Making the Hu plots
    plotting.start_plot('N964 [Mag]', 'z - N964 [Mag]')
    plt.scatter(n_mag, z_mag - n_mag, s=1, color='k', alpha=0.5)
    plt.scatter(n_mag[idx], z_mag[idx] - n_mag[idx], marker='*', s=50, color='m')
    plt.xlim(12.5, 28)
    plt.ylim(-2.4,10)
    plt.axhline(1.9, color='r', lw=1)
    plotting.end_plot('plots/hu_plot_z.png')

    plotting.start_plot('N964 [Mag]', 'i - N964 [Mag]')
    plt.scatter(n_mag, i_mag - n_mag, s=1, color='k', alpha=0.5)
    plt.scatter(n_mag[idx], i_mag[idx] - n_mag[idx], marker='*', s=50, color='m')
    plt.xlim(12.5, 28)
    plt.ylim(-2.4,10)
    plt.axhline(0.8, color='r', lw=1)
    plotting.end_plot('plots/hu_plot_i.png')
    plt.show()


    fig, ax = plt.subplots()
    #coll = ax.scatter(testData[:,0], testData[:,1], color=["blue"]*len(testData), picker = 5, s=[50]*len(testData))
    coll = ax.scatter(n_mag, z_mag - n_mag, s=1, color='k', alpha=0.5, picker = True)
    ax.scatter(n_mag[idx], z_mag[idx] - n_mag[idx], marker='*', s=100, color='m')
    ax.set_xlim(12.5, 28)
    ax.set_ylim(-2.4, 10)

    plt.axhline(1.9, color='r', lw=1)

    def on_pick(event):
        """what to do when clicking the button"""
        print('Position: ', n_ra[event.ind][0], n_dec[event.ind][0])
        ps.show_stamps(n_ra[event.ind][0], n_dec[event.ind][0])


    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()
