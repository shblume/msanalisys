""" Functions to ease plotting with mathplotlib.pyplot, statistcs, pandas and some aditional calculations. """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects

def multiplot(x_axis,y_axis, xname=None, yname=None, tname=None, note='ko', y_max=None, x_max=None,subplot=None,coord=None, y_scale='linear', a=1, y_min=None, x_min=None):
    """ Pyplot plottage for any instance. """
    # For simple plottage
    if subplot == None:
        fig, ax = plt.subplots()
        ax.plot(x_axis,y_axis,note,alpha=a)
        ax.set(xlabel=xname, ylabel=yname, title=tname)
        ax.set_yscale(y_scale)
        if y_max != None:
            plt.ylim(ymax=y_max)
        if x_max != None:
            plt.xlim(xmax=x_max)
        if y_min != None:
            plt.ylim(ymin=y_min)
        if x_min != None:
            plt.xlim(xmin=x_min)
    
    # For multiplotage
    else:
        print(' Warning: for subplots you must call "fig, ax = plt.subplots()" at least to run the code.')
        if coord == None:
            axs[subplot].plot(x_axis,y_axis,note,alpha=a)
            axs[subplot].set(xlabel=xname, ylabel=yname, title=tname)
            axs[subplot].set_yscale(y_scale)
            if y_max != None:
                plt.ylim(ymax=y_max)
            if x_max != None:
                plt.xlim(xmax=x_max)
        else:
            axs[subplot][coord].plot(x_axis,y_axis,note)
            axs[subplot][coord].set(xlabel=xname, ylabel=yname, title=tname)
            axs[subplot][coord].set_yscale(y_scale)
            if y_max != None:
                plt.ylim(ymax=y_max)
            if x_max != None:
                plt.xlim(xmax=x_max)

def center_spines(ax=None, centerx=0, centery=0,grdval=True):
    """Centers the axis spines at <centerx, centery> on the axis "ax", and
    places arrows at the end of the axis spines."""
    if ax is None:
        ax = plt.gca()

    # Set the axis's spines to be centered at the given point
    # (Setting all 4 spines so that the tick marks go in both directions)
    ax.spines['left'].set_position(('data', centerx))
    ax.spines['bottom'].set_position(('data', centery))
    ax.spines['right'].set_position(('data', centerx - 1))
    ax.spines['top'].set_position(('data', centery - 1))

    # On both the x and y axes...
    for axis, center in zip([ax.xaxis, ax.yaxis], [centerx, centery]):
        # Turn on minor and major gridlines and ticks
        axis.set_ticks_position('both')
        axis.grid(grdval, 'major', ls='solid', lw=0.5, color='gray')
        axis.grid(grdval, 'minor', ls='solid', lw=0.1, color='gray')
        axis.set_minor_locator(mpl.ticker.AutoMinorLocator())

def normalizer(n,d):
    """ Normalizes a value for another, considering the divisions by zero. Returns "None" if it goes to infinity """
    if d == 0:
        if n == 0:
            r = 0
        elif n != 0: 
            r = None 
    elif d != 0:
        r = (n/d)
    return r
