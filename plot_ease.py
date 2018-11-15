""" Functions to ease the use of the mathplotlib.pyplot library and some aditional calculations. """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects

# Functions defined:

class CenteredFormatter(matplotlib.ticker.ScalarFormatter):

    """Acts exactly like the default Scalar Formatter, but yields an empty
    label for ticks at "center"."""

    center = 0

    def __call__(self, value, pos=None):
        if value == self.center:
            return ''
        else:
            return matplotlib.ticker.ScalarFormatter.__call__(self,
value, pos)

def multiplot(x_axis,y_axis, xname=None, yname=None, tname=None, note='ko',
              y_max=None, y_min=None, x_max=None, x_min=None, alp=1,
              fig=None,axs=None,coord=None, tl=False,
              y_scale='linear',x_scale='linear'):
    """ Pyplot plottage for any instance. """
    # For simple plottage
    if fig == None:
        fig, ax = plt.subplots(tight_layout=tl)
        ax.plot(x_axis,y_axis,note,alpha=alp)
        ax.set(xlabel=xname, ylabel=yname, title=tname)
        ax.set_yscale(y_scale)
        ax.set_xscale(x_scale)
        if y_max != None:
            plt.ylim(ymax=y_max)
        if y_min != None:
            plt.ylim(ymin=y_min)
        if x_max != None:
            plt.xlim(xmax=x_max)
        if x_min != None:
            plt.xlim(xmin=x_min)
    
    # For multiplotage
    else:
        if len(coord) == 1:
            axs[coord[0]].plot(x_axis,y_axis,note,alpha=alp)
            axs[coord[0]].set(xlabel=xname, ylabel=yname, title=tname)
            axs[coord[0]].set_yscale(y_scale)
            axs[coord[0]].set_xscale(x_scale)
            if y_max != None:
                plt.ylim(ymax=y_max)
            if y_min != None:
                plt.ylim(ymin=y_min)
            if x_max != None:
                plt.xlim(xmax=x_max)
            if x_min != None:
                plt.xlim(xmin=x_min)
        else:
            axs[coord[0]][coord[1]].plot(x_axis,y_axis,note,alpha=alp)
            axs[coord[0]][coord[1]].set(xlabel=xname, ylabel=yname, title=tname)
            axs[coord[0]][coord[1]].set_yscale(y_scale)
            axs[coord[0]][coord[1]].set_xscale(x_scale)
            if y_max != None:
                plt.ylim(ymax=y_max)
            if y_min != None:
                plt.ylim(ymin=y_min)
            if x_max != None:
                plt.xlim(xmax=x_max)
            if x_min != None:
                plt.xlim(xmin=x_min)

def center_spines(ax=None, centerx=0, centery=0):
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

    # Draw an arrow at the end of the spines
    # Hide the line (but not ticks) for "extra" spines

    for side in ['right', 'top']:
        ax.spines[side].set_color('none')

    # On both the x and y axes...

    for (axis, center) in zip([ax.xaxis, ax.yaxis], [centerx, centery]):

        # Turn on minor and major gridlines and ticks

        axis.set_ticks_position('both')
        axis.grid(True, 'major', ls='solid', lw=0.5, color='gray')
        axis.grid(True, 'minor', ls='solid', lw=0.1, color='gray')
        axis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

        # Hide the ticklabels at <centerx, centery>

        formatter = CenteredFormatter()
        formatter.center = center
        axis.set_major_formatter(formatter)

    # Add offset ticklabels at <centerx, centery> using annotation
    # (Should probably make these update when the plot is redrawn...)

    (xlabel, ylabel) = map(formatter.format_data, [centerx, centery])
    ax.annotate(
        '(%s, %s)' % (xlabel, ylabel),
        (centerx, centery),
        xytext=(-4, -4),
        textcoords='offset points',
        ha='right',
        va='top',
)

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
    
""" End of the code. """
