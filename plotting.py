"""Plotting package to make all the plots look like publish ready plots."""

import matplotlib as mpl
import pylab as plt

def apply_global_settings() -> None:
    """Global settings."""
    mpl.rcParams.update({'font.size': 2})
    #mpl.rcParams['font.family'] = 'Avenir'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 2
    mpl.rc('xtick', labelsize=10)
    mpl.rc('ytick', labelsize=10)

def prettify_plot(x_label: str, y_label: str) -> plt.Figure:
    """Makes the plots look good."""
    plt.xlabel(x_label,fontsize=12)
    plt.ylabel(y_label,fontsize=12)
    plt.minorticks_on()
    plt.tick_params(which='both', width=2,direction='in')
    plt.tick_params(which='major', length=4, direction='in')

def start_plot(x_label, y_label):
    """Starting the plot."""
    fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
    prettify_plot(x_label, y_label)
    return fig


def end_plot(outfile: str) -> None:
    """Saves the figure correctly."""
    plt.savefig(outfile, bbox_inches = 'tight', transparent = False)

apply_global_settings()
