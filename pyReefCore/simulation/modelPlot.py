##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the pyReefCore synthetic coral reef core model app.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Here we set plotting functions used to visualise pyReef dataset.
"""

import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

from scipy.ndimage.filters import gaussian_filter

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

class modelPlot():
    """
    Class for plotting outputs from pyReef model.
    """

    def __init__(self, input=None):
        """
        Constructor.
        """

        self.names = input.speciesName

        return

    def speciesTime(self, pop=None, time=None, colors=None, size=(10,5), font=9, dpi=80, fname=None):
        """
        This function estimates the coral growth based on newly computed population.

        Parameters
        ----------

        variable : pop
            Population evolution for each simulated species.

        variable : time
            Time step evolution

        variable : colors
            Matplotlib color map to use

        variable : size
            Figure size

        variable : font
            Figure font size

        variable : dpi
            Figure resolution

        variable : fname
            Save PNG filename.
        """

        matplotlib.rcParams.update({'font.size': font})

        # Define figure size
        fig, ax = plt.subplots(1,figsize=size, dpi=dpi)
        ax.set_axis_bgcolor('#f2f2f3')

        # Plotting curves
        for s in range(len(pop)):
            ax.plot(time, pop[s,:], label=self.names[s],linewidth=3,c=colors[s])

        # Legend, title and labels
        plt.grid()
        lgd = plt.legend(frameon=False,loc=4,prop={'size':font+1}, bbox_to_anchor=(1.2,-0.02))
        plt.xlabel('Time [y]',size=font+2)
        plt.ylabel('Population',size=font+2)
        plt.ylim(0., int(pop.max())+1)


        ttl = ax.title
        ttl.set_position([.5, 1.05])
        plt.title('Evolution of species populations with time',size=font+3)
        plt.show()

        if fname is not None:
            fig.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')

        return

    def speciesDepth(self, pop=None, depth=None, colors=None, size=(10,5), font=9, dpi=80, fname=None):
        """
        Variation of coral growth with depth

        Parameters
        ----------

        variable : pop
            Population evolution for each simulated species.

        variable : depth
            Depth of coral core

        variable : colors
            Matplotlib color map to use

        variable : size
            Figure size

        variable : font
            Figure font size

        variable : dpi
            Figure resolution

        variable : fname
            Save PNG filename.
        """

        matplotlib.rcParams.update({'font.size': font})

        # Define figure size
        fig, ax = plt.subplots(1,figsize=size, dpi=dpi)
        ax.set_axis_bgcolor('#f2f2f3')

        # Plotting curves
        for s in range(len(pop)):
            ax.plot(np.cumsum(depth), pop[s,1:], label=self.names[s],linewidth=3,c=colors[s])

        # Legend, title and labels
        plt.grid()
        lgd = plt.legend(frameon=False,loc=4,prop={'size':font+1}, bbox_to_anchor=(1.2,-0.02))
        plt.xlabel('Depth [m]',size=font+2)
        plt.ylabel('Population',size=font+2)
        plt.ylim(0., int(pop.max())+1)


        ttl = ax.title
        ttl.set_position([.5, 1.05])
        plt.title('Evolution of species populations with depth',size=font+3)
        plt.show()

        if fname is not None:
            fig.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')

        return

    def drawCore(self, pop=None, depth=None, depthext = None, thext = None, propext = [0.,1.], surf=None, lwidth = 3,
                 colors=None, size=(8,10), font=8, dpi=80, figname=None, filename = None, sep = '\t'):
        """
        Plot core evolution

        Parameters
        ----------

        variable : pop
            Population evolution for each simulated species.

        variable : depth
            Depth of coral core

        variable : surf
            Surface elevation of the coral system

        variable : depthext
            Core depth extension to plot [m]

        variable : thext
            Core thickness range to plot [m]

        variable : propext
            Core ranging proportion to plot between [0,1.]

        variable : lwidth
            Figure lines width

        variable : colors
            Matplotlib color map to use

        variable : size
            Figure size

        variable : font
            Figure font size

        variable : dpi
            Figure resolution

        variable : figname
            Save gigure (the type of file needs to be provided e.g. .png or .pdf).

        variable : filename
            Save model output to a CSV file.

        variable : sep
            Separator used in the CSV file.
        """

        p1 = pop
        p2 = pop/depth
        p3 = np.cumsum(pop/depth,axis=0)
        d = 2.*surf - np.cumsum(depth)

        # Define figure size
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=size, sharey=True, dpi=dpi)
        ax3.set_axis_bgcolor('#f2f2f3')

        # Plotting curves
        for s in range(len(pop)):
            ax1.plot(p1[s,:], d, label=self.names[s], linewidth=lwidth, c=colors[s])
            ax2.plot(p2[s,:], d, label=self.names[s], linewidth=lwidth, c=colors[s])
            if s == 0:
                ax3.fill_betweenx(d, 0, p3[s,:], color=colors[s])
            else:
                ax3.fill_betweenx(d, p3[s-1,:], p3[s,:], color=colors[s])
            ax3.plot(p3[s,:], d, 'k--', label=self.names[s], linewidth=lwidth-1)

        # Legend, title and labels
        ax1.grid()
        ax2.grid()
        ax3.grid()
        lgd = ax1.legend(frameon=False, loc=1, prop={'size':font+1}, bbox_to_anchor=(4.,0.2))
        ax1.locator_params(axis='x', nbins=4)
        ax2.locator_params(axis='x', nbins=5)
        ax3.locator_params(axis='x', nbins=5)
        ax1.locator_params(axis='y', nbins=10)

        # Axis
        ax1.set_ylabel('Depth below ocean surface [m]', size=font+4)
        ax1.set_ylim(depthext[1], depthext[0])
        ax1.set_xlim(thext[0], thext[1])
        ax2.set_ylim(depthext[1], depthext[0])
        ax2.set_xlim(propext[0], propext[1])
        ax3.set_ylim(depthext[1], depthext[0])
        ax3.set_xlim(0., 1.)
        ax1.xaxis.tick_top()
        ax2.xaxis.tick_top()
        ax3.xaxis.tick_top()
        ax1.tick_params(axis='y', pad=5)
        ax1.tick_params(axis='x', pad=5)
        ax2.tick_params(axis='x', pad=5)
        ax3.tick_params(axis='x', pad=5)

        # Title
        tt1 = ax1.set_title('Thickness [m]', size=font+3)
        tt2 = ax2.set_title('Proportion [%]', size=font+3)
        tt3 = ax3.set_title('Accumulated [%]', size=font+3)
        tt1.set_position([.5, 1.04])
        tt2.set_position([.5, 1.04])
        tt3.set_position([.5, 1.04])
        fig.tight_layout()

        plt.show()

        if figname is not None:
            fig.savefig(figname, bbox_extra_artists=(lgd,), bbox_inches='tight')

        if filename is not None:
            tmp = np.column_stack((d.T,p1.T))
            tmp1 = np.column_stack((tmp,p2.T))
            tmp2 = np.column_stack((tmp1,p3.T))

            cols = []
            cols.append('depth')
            for s in range(len(self.names)):
                cols.append('th_'+self.names[s])
            for s in range(len(self.names)):
                cols.append('prop_'+self.names[s])
            for s in range(len(self.names)):
                cols.append('acc_'+self.names[s])

            df = pd.DataFrame(tmp2)
            df.columns = cols
            df.to_csv(filename, sep=sep, encoding='utf-8', index=False)
            print 'Model results have been saved in',filename

        return
