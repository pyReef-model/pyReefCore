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

        self.names = np.empty(input.speciesNb+1, dtype="S14")
        self.names[:input.speciesNb] = input.speciesName
        self.names[-1] = 'silicilastic'
        self.step = int(input.laytime/input.tCarb)
        self.pop = None
        self.timeCarb = None
        self.depth = None
        self.accspace = None
        self.timeLay = None
        self.surf = None
        self.sedH = None
        self.sealevel = None
        self.mbsl = None
        self.sedinput = None
        self.waterflow = None
        self.folder = input.outDir

        return

    def two_scales(self, ax1, time, data0, data1, c1, c2, font):
        """

        Parameters
        ----------
        ax : axis
            Axis to put two scales on

        time : array-like
            x-axis values for both datasets

        data1: array-like
            Data for left hand scale

        data2 : array-like
            Data for right hand scale

        c1 : color
            Color for line 1

        c2 : color
            Color for line 2

        Returns
        -------
        ax : axis
            Original axis
        ax2 : axis
            New twin axis
        """
        ax2 = ax1.twinx()

        ax1.plot(time, data0, color=c1, linewidth=3)
        ax1.set_xlabel('Time [y]',size=font+2)
        ax1.set_ylabel('accomodation space [m]',size=font+2)
        ax1.yaxis.label.set_color(c1)

        ax2.plot(time, data1, color=c2, linewidth=3)
        ax2.set_ylabel('water depth [mbsl]',size=font+2)
        ax2.yaxis.label.set_color(c2)

        return ax1, ax2

    def color_y_axis(self, ax, color):
        """Color your axes."""
        for t in ax.get_yticklabels():
            t.set_color(color)

        return None

    def accomodationTime(self, colors=None, size=(10,5), font=9, dpi=80, fname=None):
        """
        This function estimates the accomodation space through time.

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

        if colors is not None:
            c1 = colors[0]
            c2 = colors[-1]
        else:
            c1 = '#1f77b4'
            c2 = '#ff7f0e'

        # Define figure size
        fig, ax = plt.subplots(1,figsize=size, dpi=dpi)
        ax.set_facecolor('#f2f2f3')
        tmp = self.mbsl[:-2]-self.accspace[:-2]


        tmp2 = np.ediff1d(tmp)
        d = np.zeros(len(self.accspace[:-2]))
        d[1:] = tmp2
        d[0] = d[1]
        d[-1] = d[-2]

        ax1, ax2 = self.two_scales(ax,self.timeCarb[:-2],self.accspace[:-2],tmp,c1,c2,font)

        # Plotting curves
        #ax.plot(self.timeCarb[:-2], self.accspace[:-2], linewidth=3,c=colors)
        #ax.plot(self.timeCarb[:-2], tmp, linewidth=3,c=colors)

        #plt.xlabel('Time [y]',size=font+2)
        #plt.ylabel('accomodation space [m]',size=font+2)

        ttl = ax.title
        ttl.set_position([.5, 1.05])
        plt.title('Accomodation space & water depth evolution through time',size=font+3)

        self.color_y_axis(ax1, c1)
        self.color_y_axis(ax2, c2)

        plt.xlim(0., self.timeCarb.max())

        # Legend, title and labels
        plt.grid()
        plt.show()

        if fname is not None:
            name = self.folder+'/'+fname
            fig.savefig(name, bbox_inches='tight')

        # Define figure size
        fig, ax = plt.subplots(1,figsize=size, dpi=dpi)
        ax.set_facecolor('#f2f2f3')

        # Plotting curves
        ax.plot(self.timeCarb[:-2], d, linewidth=3,c='#2ca02c')

        ttl = ax.title
        ttl.set_position([.5, 1.05])
        plt.title('Core production thickness through time',size=font+3)

        plt.xlabel('Time [y]',size=font+2)
        plt.ylabel('Core thickness [m]',size=font+2)
        ax.yaxis.label.set_color('#2ca02c')
        plt.xlim(0., self.timeCarb.max())

        # Legend, title and labels
        plt.grid()
        plt.show()

        if fname is not None:
            name = self.folder+'/prodvsdepth-'+fname
            fig.savefig(name, bbox_inches='tight')

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
        ax.set_facecolor('#f2f2f3')

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
            name = self.folder+'/'+fname
            fig.savefig(name, bbox_extra_artists=(lgd,), bbox_inches='tight')

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
        ax.set_facecolor('#f2f2f3')

        # Plotting curves
        bottom = self.surf + self.depth.sum()
        d = bottom - np.cumsum(self.depth)
        for s in range(len(self.pop)):
            ax.plot(d, self.pop[s,::self.step], label=self.names[s],linewidth=3,c=colors[s])

        # Legend, title and labels
        plt.grid()
        lgd = plt.legend(frameon=False,loc=4,prop={'size':font+1}, bbox_to_anchor=(1.2,-0.02))
        plt.xlabel('Depth [m]',size=font+2)
        plt.ylabel('Population',size=font+2)
        plt.ylim(0., int(self.pop.max())+1)
        plt.xlim(d.max(), d.min())

        ttl = ax.title
        ttl.set_position([.5, 1.05])
        plt.title('Evolution of species populations with depth',size=font+3)
        plt.show()

        if fname is not None:
            name = self.folder+'/'+fname
            fig.savefig(name, bbox_extra_artists=(lgd,), bbox_inches='tight')

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

        p1 = self.sedH[:,:-1]
        ids = np.where(self.depth[:-1]>0)[0]
        p2 = np.zeros((self.sedH.shape))
        p3 = np.zeros((self.sedH.shape))

        p2[:,ids] = self.sedH[:,ids]/self.depth[ids]
        p3[:,ids] = np.cumsum(self.sedH[:,ids]/self.depth[ids],axis=0)
        bottom = self.surf + self.depth[:-1].sum()
        d = bottom - np.cumsum(self.depth[:-1])
        facies = np.argmax(p1, axis=0)

        if thext == None:
            thext = [0.,p1.max()]

        if depthext == None:
            depthext = [self.surf,bottom-self.depth[0]]


        colsed[len(self.sedH)-1]=np.array([244./256.,164/256.,96/256.,1.])

        # Define figure size
        fig = plt.figure(figsize=size, dpi=dpi)
        gs = gridspec.GridSpec(1,19)
        ax1 = fig.add_subplot(gs[:5])
        ax2 = fig.add_subplot(gs[5:10], sharey=ax1)
        ax3 = fig.add_subplot(gs[10:15], sharey=ax1)
        ax42 = fig.add_subplot(gs[15:17], frame_on=False)
        ax4 = fig.add_subplot(gs[15:17], sharey=ax1)
        ax52 = fig.add_subplot(gs[17:19], frame_on=False)
        ax5 = fig.add_subplot(gs[17:19], sharey=ax1)
        ax3.set_facecolor('#f2f2f3')
        ax4.set_facecolor('#f2f2f3')
        ax5.set_facecolor('#f2f2f3')
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
            ax2.plot(p2[s,:-1], d, label=self.names[s], linewidth=lwidth, c=colsed[s])
            if s == 0:
                ax3.fill_betweenx(d, 0, p3[s,:-1], color=colsed[s])
            else:
                ax3.fill_betweenx(d, p3[s-1,:-1], p3[s,:-1], color=colsed[s])
            ax3.plot(p3[s,:-1], d, 'k--', label=self.names[s], linewidth=lwidth-1)

        tmpx = np.zeros(len(self.timeLay))
        tmpx[-1] = 1
        ax42.plot(tmpx[:-1], self.timeLay[:-1], zorder=1)
        ax42.yaxis.tick_right()
        ax52.plot(tmpx[:-1], self.timeLay[:-1], zorder=1)
        ax52.yaxis.tick_right()

        for s in range(len(d)):
            y[0] = d[s]
            y[1] = d[s]
            ax4.fill_between(x, old, y, color=coltime[s], zorder=10)
            ax5.fill_between(x, old, y, color=colsed[facies[s]], zorder=10)
            old[0] = y[0]
            old[1] = y[1]
            ax4.plot(x,y,'k', zorder=10,linewidth=1)
            ax5.plot(x,y,'k', zorder=10,linewidth=1)

        # Legend, title and labels
        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax42.get_xaxis().set_visible(False)
        ax4.get_xaxis().set_visible(False)
        #ax4.get_yaxis().set_visible(False)
        ax52.get_xaxis().set_visible(False)
        ax5.get_xaxis().set_visible(False)
        #ax5.get_yaxis().set_visible(False)
        # ax5.get_yaxis().set_visible(False)
        lgd = ax1.legend(frameon=False, loc=1, prop={'size':font+1}, bbox_to_anchor=(6.2,0.2))
        ax1.locator_params(axis='x', nbins=5)
        ax2.locator_params(axis='x', nbins=5)
        ax3.locator_params(axis='x', nbins=5)
        ax1.locator_params(axis='y', nbins=10)

        # Axis
        ax1.set_ylabel('Depth below present mean sea-level [m]', size=font+4)
        ax1.set_ylim(depthext[1], depthext[0])
        ax1.set_xlim(thext[0], thext[1]+thext[1]*0.1)
        ax2.set_ylim(depthext[1], depthext[0])
        ax2.set_xlim(propext[0], propext[1])
        ax3.set_ylim(depthext[1], depthext[0])
        ax42.set_ylim(self.timeLay[0],self.timeLay[-1])
        ax4.set_ylim(depthext[1], depthext[0])
        ax52.set_ylim(self.timeLay[0],self.timeLay[-1])
        ax5.set_ylim(depthext[1], depthext[0])
        ax3.set_xlim(0., 1.)
        ax1.xaxis.tick_top()
        ax2.xaxis.tick_top()
        ax3.xaxis.tick_top()
        ax1.tick_params(axis='y', pad=5)
        #ax4.tick_params(axis='y', pad=5)
        ax42.tick_params(axis='y', pad=5)
        #ax42.set_ylabel('Time [years]', size=font+4)
        #ax42.yaxis.set_label_position('right')

        ax52.tick_params(axis='y', pad=5)
        ax5.tick_params(axis='y', pad=5)
        #ax5.set_ylabel('Time [years]', size=font+4)
        ax5.yaxis.set_label_position('right')

        ax1.tick_params(axis='x', pad=5)
        ax2.tick_params(axis='x', pad=5)
        ax3.tick_params(axis='x', pad=5)

        # Title
        tt1 = ax1.set_title('Thickness [m]', size=font+3)
        tt2 = ax2.set_title('Proportion [%]', size=font+3)
        tt3 = ax3.set_title('Accumulated [%]', size=font+3)
        tt4 = ax4.set_title('Time \n layers', size=font+3)
        tt5 = ax5.set_title('Bio. \n facies', size=font+3)
        tt1.set_position([.5, 1.04])
        tt2.set_position([.5, 1.04])
        tt3.set_position([.5, 1.04])
        tt4.set_position([.5, 1.025])
        tt5.set_position([.5, 1.025])
        fig.tight_layout()
        plt.tight_layout()
        #labels = [item.get_text() for item in ax2.get_yticklabels()]
        #for l in range(len(labels)):
        #    labels[l] = ' '
        #ax2.set_yticklabels(labels)
        #ax3.set_yticklabels(labels)
        plt.figtext(0.885, 0.005, 'left axis:depth [m] - right axis:time [years]',horizontalalignment='center', fontsize=font)
        plt.show()


        if figname is not None:
            name = self.folder+'/'+figname
            fig.savefig(name, bbox_extra_artists=(lgd,), bbox_inches='tight')
            print 'Figure has been saved in',name

        # Define figure size
        fig = plt.figure(figsize=size, dpi=dpi)
        gs = gridspec.GridSpec(1,11)
        ax1 = fig.add_subplot(gs[:3])
        ax2 = fig.add_subplot(gs[3:6], sharey=ax1)
        ax3 = fig.add_subplot(gs[6:9], sharey=ax1)

        ax1.plot(self.sealevel, self.timeLay, linewidth=lwidth, c='slateblue')
        ax2.plot(self.waterflow, self.timeLay, linewidth=lwidth, c='darkcyan')
        ax3.plot(self.sedinput, self.timeLay, linewidth=lwidth, c='sandybrown')

        ax1.set_ylabel('Simulation time [a]', size=font+4)
        ax1.set_ylim(self.timeLay.min(), self.timeLay.max())
        ax2.set_ylim(self.timeLay.min(), self.timeLay.max())
        ax3.set_ylim(self.timeLay.min(), self.timeLay.max())
        ax1.set_facecolor('#f2f2f3')
        ax2.set_facecolor('#f2f2f3')
        ax3.set_facecolor('#f2f2f3')
        ax1.locator_params(axis='x', nbins=5)
        ax2.locator_params(axis='x', nbins=5)
        ax3.locator_params(axis='x', nbins=5)

        # Title
        tt1 = ax1.set_title('Sea-level [m]', size=font+3)
        tt2 = ax2.set_title('Water flow [m/d]', size=font+3)
        tt3 = ax3.set_title('Sediment input [m/d]', size=font+3)
        tt1.set_position([.5, 1.04])
        tt2.set_position([.5, 1.04])
        tt3.set_position([.5, 1.04])
        fig.tight_layout()
        plt.tight_layout()
        plt.show()
        if figname is not None:
            name = self.folder+'/envi'+figname
            fig.savefig(name, bbox_extra_artists=(lgd,), bbox_inches='tight')
            print 'Figure has been saved in','envi'+name
        print ''

        if filename is not None:
            tmp = np.column_stack((d.T,p1.T))
            tmp1 = np.column_stack((tmp,p2[:,:-1].T))
            tmp2 = np.column_stack((tmp1,p3[:,:-1].T))
            tmp3 = np.column_stack((tmp2,self.sealevel[:-1].T))
            tmp4 = np.column_stack((tmp3,self.waterflow[:-1].T))
            tmp5 = np.column_stack((tmp4,self.sedinput[:-1].T))

            cols = []
            cols.append('depth')
            for s in range(len(self.names)):
                cols.append('th_'+self.names[s])
            for s in range(len(self.names)):
                cols.append('prop_'+self.names[s])
            for s in range(len(self.names)):
                cols.append('acc_'+self.names[s])
            cols.append('sealevel')
            cols.append('waterflow')
            cols.append('sedinput')

            df = pd.DataFrame(tmp5)
            df.columns = cols
            name = self.folder+'/'+filename
            df.to_csv(name, sep=sep, encoding='utf-8', index=False)
            print 'Model results have been saved in',name

        return
