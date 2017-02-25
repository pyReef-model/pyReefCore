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

        # Core parameters size based on layer number
        self.layNb = int((input.tEnd - input.tStart)/input.laytime)
        self.thickness = numpy.zeros(self.layNb,dtype=float)
        self.coralH = numpy.zeros((input.speciesNb,self.layNb),dtype=float)

        # Diagonal part of the community matrix (coefficient ii)
        self.alpha = input.communityMatrix.diagonal()
        self.layTime = numpy.arange(input.tStart, input.tEnd+input.laytime, input.laytime)

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
        # print self.prod
        # print self.alpha
        # print epsilon
        # print coral
        # print self.dt
        # print toth
        # print production
        # print layID
        # print ''

        return
