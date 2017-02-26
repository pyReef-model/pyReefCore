##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the pyReefCore synthetic coral reef core model app.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module builds the core records through time based on coral species evolution and
the interactions between the active forcing paramters.
"""
import os
import numpy
import pandas as pd
import skfuzzy as fuzz

import matplotlib
import matplotlib.pyplot as plt

class coreData:
    """
    This class defines the core parameters
    """

    def __init__(self, input = input):
        """
        Constructor.
        """

        self.dt = input.tCarb

        # Initial core depth
        self.topH = input.depth0

        # Production rate for each carbonate
        self.prod = input.speciesProduction
        self.names = input.speciesName

        # Core parameters size based on layer number
        self.layNb = int((input.tEnd - input.tStart)/input.laytime)
        self.thickness = numpy.zeros(self.layNb,dtype=float)
        self.coralH = numpy.zeros((input.speciesNb,self.layNb),dtype=float)

        # Diagonal part of the community matrix (coefficient ii)
        self.communityMatrix = input.communityMatrix
        self.alpha = input.communityMatrix.diagonal()
        self.layTime = numpy.arange(input.tStart, input.tEnd+input.laytime, input.laytime)

        # Shape functions
        self.edepth = input.enviDepth
        self.eflow = input.enviFlow
        self.esed = input.enviSed

        return

    def _plot_fuzzy_curve(self, xd, xs, xf, dtrap, strap, ftrap, size,
                          dpi, font, colors, width, fname):

        matplotlib.rcParams.update({'font.size': font})

        for s in range(len(self.names)):
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=size, sharey=True, dpi=dpi)
            ax1.set_axis_bgcolor('#f2f2f3')
            ax2.set_axis_bgcolor('#f2f2f3')
            ax3.set_axis_bgcolor('#f2f2f3')
            fig.tight_layout()
            ax1.grid()
            ax2.grid()
            ax3.grid()
            ax1.plot(xd, dtrap[s], linewidth=width, label=self.names[s],c=colors[s])
            ax2.plot(xs, strap[s], linewidth=width, label=self.names[s],c=colors[s])
            ax3.plot(xf, ftrap[s], linewidth=width, label=self.names[s],c=colors[s])
            ax1.set_ylabel('Environmental Factor',size=font+3)
            ax1.set_ylim(-0.1, 1.1)
            ax1.set_xlabel('Water Depth [m]',size=font+2)
            ax3.set_xlabel('Flow Velocity [m/d]',size=font+2)
            ax3.set_ylim(-0.1, 1.1)
            ax2.set_xlabel('Sediment Input [m/d]',size=font+2)
            ax2.set_ylim(-0.1, 1.1)
            ax3.yaxis.set_label_position("right")
            ax3.set_ylabel(self.names[s],size=font+3,fontweight='bold')
            plt.show()
            if fname is not None:
                fig.savefig(name[s]+fname)

        return

    def initialSetting(self, font=8, size=(8,2.5), width=3, dpi=80, fname=None):
        """
        Visualise the initial conditions of the model run.

        Parameters
        ----------
        variable : font
            Environmental shape figures font size

        variable : size
            Environmental shape figures size

        variable : width
            Environmental shape figures line width

        variable : dpi
            Figure resolution

        variable : fname
            Save filename.
        """

        from matplotlib.cm import terrain
        nbcolors = len(self.names)+3
        colors = terrain(numpy.linspace(0, 1, nbcolors))

        print 'Community matrix aij representing the interactions between species:'
        print ''
        cols = []
        ids = []
        for i in range(len(self.names)):
            cols.append('a'+str(i))
            ids.append('a'+str(i)+'j')
        df = pd.DataFrame(self.communityMatrix, index=ids)
        df.columns = cols
        print df
        print ''
        print 'Species maximum production rates [m/y]:'
        print ''
        index = [self.names]
        df = pd.DataFrame(self.prod,index=index)
        df.columns = ['Prod.']
        print df
        print ''
        print 'Environmental trapezoidal shape functions:'

        # Visualise fuzzy production curve
        xs = numpy.linspace(0, self.esed.max(), num=201, endpoint=True)
        xf = numpy.linspace(0, self.eflow.max(), num=201, endpoint=True)
        xd = numpy.linspace(0, self.edepth.max(), num=201, endpoint=True)
        dtrap = []
        strap = []
        ftrap = []
        for s in range(0,len(self.names)):
            dtrap.append(fuzz.trapmf(xd, self.edepth[s,:]))
            strap.append(fuzz.trapmf(xs, self.esed[s,:]))
            ftrap.append(fuzz.trapmf(xf, self.eflow[s,:]))
        self._plot_fuzzy_curve(xd, xs, xf, dtrap, strap, ftrap, size,
                               dpi, font, colors, width, fname)

        return

    def coralProduction(self, layID, coral, epsilon):
        """
        This function estimates the coral growth based on newly computed population.

        Parameters
        ----------

        variable : layID
            Index of current stratigraphic layer.

        variable : coral
            Species population distribution at current time step.

        variable : epsilon
            Intrinsic rate of a population species (malthus parameter)
        """

        # Compute production for the given time step [m]
        production = - self.prod * self.alpha / epsilon * coral * self.dt

        # Total thickness deposited
        toth = production.sum()

        # Update current layer composition
        self.coralH[:,layID] += production
        # Update current layer thickness
        self.thickness[layID] += toth
        # Update current layer top elevation
        self.topH -= toth

        return
