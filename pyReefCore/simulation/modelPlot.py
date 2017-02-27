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
from matplotlib import gridspec
import matplotlib.pyplot as plt

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
        self.step = int(input.laytime/input.tCarb)
        self.pop = None
        self.timeCarb = None
        self.depth = None
        self.timeLay = None
        self.surf = None
        self.sedH = None

        return

    def speciesTime(self, colors=None, size=(10,5), font=9, dpi=80, fname=None):
        """
        This function estimates the coral growth based on newly computed population.

        Parameters
        ----------

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
        for s in range(len(self.pop)):
            ax.plot(self.timeCarb, self.pop[s,:], label=self.names[s],linewidth=3,c=colors[s])

        # Legend, title and labels
        plt.grid()
        lgd = plt.legend(frameon=False,loc=4,prop={'size':font+1}, bbox_to_anchor=(1.2,-0.02))
        plt.xlabel('Time [y]',size=font+2)
        plt.ylabel('Population',size=font+2)
        plt.ylim(0., int(self.pop.max())+1)
        plt.xlim(0., self.timeCarb.max())


        ttl = ax.title
        ttl.set_position([.5, 1.05])
        plt.title('Evolution of species populations with time',size=font+3)
        plt.show()

        if fname is not None:
            fig.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')

        return

    def speciesDepth(self, colors=None, size=(10,5), font=9, dpi=80, fname=None):
        """
        Variation of coral growth with depth

        Parameters
        ----------

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
        for s in range(len(self.pop)):
            ax.plot(np.cumsum(self.depth), self.pop[s,::self.step], label=self.names[s],linewidth=3,c=colors[s])

        # Legend, title and labels
        plt.grid()
        lgd = plt.legend(frameon=False,loc=4,prop={'size':font+1}, bbox_to_anchor=(1.2,-0.02))
        plt.xlabel('Depth [m]',size=font+2)
        plt.ylabel('Population',size=font+2)
        plt.ylim(0., int(self.pop.max())+1)
        plt.xlim(0., self.depth.sum())

        ttl = ax.title
        ttl.set_position([.5, 1.05])
        plt.title('Evolution of species populations with depth',size=font+3)
        plt.show()

        if fname is not None:
            fig.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')

        return

    def drawCore(self, depthext = None, thext = None, propext = [0.,1.], lwidth = 3,
                 colsed=None, coltime=None, size=(8,10), font=8, dpi=80, figname=None,
                 filename = None, sep = '\t'):
        """
        Plot core evolution

        Parameters
        ----------

        variable : depthext
            Core depth extension to plot [m]

        variable : thext
            Core thickness range to plot [m]

        variable : propext
            Core ranging proportion to plot between [0,1.]

        variable : lwidth
            Figure lines width

        variable : colsed
            Matplotlib color map to use for production plots

        variable : coltime
            Matplotlib color map to use for time layer plots

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

        p1 = self.sedH
        p2 = self.sedH/self.depth
        p3 = np.cumsum(self.sedH/self.depth,axis=0)
        bottom = self.surf + self.depth.sum()
        d = bottom - np.cumsum(self.depth)

        if thext == None:
            thext = [0.,p1.max()]

        if depthext == None:
            depthext = [self.surf,bottom]

        # Define figure size
        fig = plt.figure(figsize=size, dpi=dpi)
        gs = gridspec.GridSpec(1,11)
        ax1 = fig.add_subplot(gs[:3])
        ax2 = fig.add_subplot(gs[3:6], sharey=ax1)
        ax3 = fig.add_subplot(gs[6:9], sharey=ax1)
        ax4 = fig.add_subplot(gs[9], sharey=ax1)
        ax5 = fig.add_subplot(gs[10], sharey=ax1)
        ax3.set_axis_bgcolor('#f2f2f3')
        ax4.set_axis_bgcolor('#f2f2f3')
        ax5.set_axis_bgcolor('#f2f2f3')
        x = np.zeros(2)
        y = np.zeros(2)
        x[0] = 0.
        x[1] = 1.
        old = np.zeros(2)
        old[0] = bottom
        old[1] = bottom

        # Plotting curves
        for s in range(len(self.sedH)):
            ax1.plot(p1[s,:], d, label=self.names[s], linewidth=lwidth, c=colsed[s])
            ax2.plot(p2[s,:], d, label=self.names[s], linewidth=lwidth, c=colsed[s])
            if s == 0:
                ax3.fill_betweenx(d, 0, p3[s,:], color=colsed[s])
            else:
                ax3.fill_betweenx(d, p3[s-1,:], p3[s,:], color=colsed[s])
            ax3.plot(p3[s,:], d, 'k--', label=self.names[s], linewidth=lwidth-1)

        for s in range(len(d)):
            y[0] = d[s]
            y[1] = d[s]
            ax4.fill_between(x, old, y, color=coltime[s])
            old[0] = y[0]
            old[1] = y[1]
            ax4.plot(x,y,'w')

        # Legend, title and labels
        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax4.get_xaxis().set_visible(False)
        ax4.get_yaxis().set_visible(False)
        ax5.get_xaxis().set_visible(False)
        ax5.get_yaxis().set_visible(False)
        lgd = ax1.legend(frameon=False, loc=1, prop={'size':font+1}, bbox_to_anchor=(5.,0.2))
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
        ax4.set_ylim(depthext[1], depthext[0])
        ax5.set_ylim(depthext[1], depthext[0])
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
        tt4 = ax4.set_title('Time \n layers', size=font+3)
        tt5 = ax5.set_title('Auto. \n facies', size=font+3)
        tt1.set_position([.5, 1.04])
        tt2.set_position([.5, 1.04])
        tt3.set_position([.5, 1.04])
        tt4.set_position([.5, 1.025])
        tt5.set_position([.5, 1.025])
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
