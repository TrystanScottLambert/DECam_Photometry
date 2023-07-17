"""Plotting package to make all the plots look like publish ready plots."""

import matplotlib as mpl
import pylab as plt

def apply_global_settings() -> None:
    """Global settings."""
    mpl.rcParams.update({'font.size': 2})
    #mpl.rcParams['font.family'] = 'Avenir'
    plt.rcParams['font.size'] = 8
    plt.rcParams['axes.linewidth'] = 1.1
    mpl.rc('xtick', labelsize=8)
    mpl.rc('ytick', labelsize=8)

def prettify_plot(x_label: str, y_label: str) -> None:
    """Makes the plots look good."""
    plt.xlabel(x_label,fontsize=12)
    plt.ylabel(y_label,fontsize=12)
    plt.minorticks_on()
    plt.tick_params(which='both', width=1.2,direction='in')
    plt.tick_params(which='major', length=3, direction='in')

def start_plot(x_label: str, y_label: str) -> plt.Figure:
    """Starting the plot."""
    fig = plt.figure(figsize = (3.54, 3.54), dpi = 600)
    prettify_plot(x_label, y_label)
    return fig

def prettify_axis(axis: plt.axes, x_label: str, y_label: str) -> None:
    """Does prettify plot but for a given axis"""
    axis.set_xlabel(x_label, fontsize=12)
    axis.set_ylabel(y_label, fontsize=12)
    axis.minorticks_on()
    axis.tick_params(which='both', width=1.2,direction='in')
    axis.tick_params(which='major', length=3, direction='in')


def end_plot(outfile: str) -> None:
    """Saves the figure correctly."""
    plt.savefig(outfile, bbox_inches = 'tight', transparent = False)

apply_global_settings()
