##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the pyReefCore synthetic coral reef core model app.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines several functions used to force pyReefCore simulation with external
processes related to sediment input, flow velocity and sea level.
"""

import os
import numpy
import pandas
import skfuzzy as fuzz
from scipy import interpolate

class enviForce:
    """
    This class defines external forcing parameters.
    """

    def __init__(self, input):
        """
        Constructor.

        Parameters
        ----------
        class: input
            Input parameter class.
        """

        self.sea0 = input.seaval
        self.seafile = input.seafile
        self.sealevel = None
        self.seatime = None
        self.seaFunc = None

        self.sed0 = input.sedval
        self.sedfile = input.sedfile
        self.sedlevel = None
        self.sedtime = None
        self.sedFunc = None

        self.flow0 = input.flowval
        self.flowfile = input.flowfile
        self.flowlevel = None
        self.flowtime = None
        self.flowFunc = None

        if self.seafile != None:
            self._build_Sea_function()

        if self.sedfile != None:
            self._build_Sed_function()

        if self.flowfile != None:
            self._build_Flow_function()

        # Shape functions
        self.edepth = None
        self.xd = None
        self.dtrap = []
        if input.seaOn:
            self.edepth = input.enviDepth
            # Trapeizoidal environment depth production curve
            self.xd = numpy.linspace(0, self.edepth.max(), num=1001, endpoint=True)
            for s in range(input.speciesNb):
                self.dtrap.append(fuzz.trapmf(self.xd, self.edepth[s,:]))

        self.speciesNb = input.speciesNb
        self.eflow = None
        self.xf = None
        self.ftrap = []
        if input.flowOn:
            self.eflow = input.enviFlow
            # Trapeizoidal environment flow production curve
            self.xf = numpy.linspace(0, self.eflow.max(), num=1001, endpoint=True)
            for s in range(input.speciesNb):
                self.ftrap.append(fuzz.trapmf(self.xf, self.eflow[s,:]))

        self.esed = None
        self.xs = None
        self.strap = []
        if input.sedOn:
            self.esed = input.enviSed
            # Trapeizoidal environment sediment production curve
            self.xs = numpy.linspace(0, self.esed.max(), num=1001, endpoint=True)
            for s in range(input.speciesNb):
                self.strap.append(fuzz.trapmf(self.xs, self.esed[s,:]))

        return

    def _extract_enviParam(self, x, xmf, xx):
        """
        Find the degree of membership ``u(xx)`` for a given value of ``x = xx``.
        """

        # Nearest discrete x-values
        x1 = x[x <= xx][-1]
        x2 = x[x >= xx][0]

        idx1 = numpy.nonzero(x == x1)[0][0]
        idx2 = numpy.nonzero(x == x2)[0][0]

        xmf1 = xmf[idx1]
        xmf2 = xmf[idx2]

        if x1 == x2:
            xxmf = xmf[idx1]
        else:
            slope = (xmf2 - xmf1) / float(x2 - x1)
            xxmf = slope * (xx - x1) + xmf1

        return xxmf

    def _build_Sea_function(self):
        """
        Using Pandas library to read the sea level file and define sea level interpolation
        function based on Scipy 1D cubic function.
        """

        # Read sea level file
        seadata = pandas.read_csv(self.seafile, sep=r'\s+', engine='c',
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)

        self.seatime = seadata.values[:,0]
        tmp = seadata.values[:,1]
        self.seaFunc = interpolate.interp1d(self.seatime, tmp, kind='cubic')

        return

    def _build_Sed_function(self):
        """
        Using Pandas library to read the sediment input file and define interpolation
        function based on Scipy 1D cubic function.
        """

        # Read sea level file
        seddata = pandas.read_csv(self.sedfile, sep=r'\s+', engine='c',
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)

        self.sedtime = seddata.values[:,0]
        tmp = seddata.values[:,1]
        self.sedFunc = interpolate.interp1d(self.sedtime, tmp, kind='cubic')

        return

    def _build_Flow_function(self):
        """
        Using Pandas library to read the flow velocity file and define interpolation
        function based on Scipy 1D cubic function.
        """

        # Read sea level file
        flowdata = pandas.read_csv(self.flowfile, sep=r'\s+', engine='c',
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)

        self.flowtime = flowdata.values[:,0]
        tmp = flowdata.values[:,1]
        self.flowFunc = interpolate.interp1d(self.flowtime, tmp, kind='cubic')

        return

    def getSea(self, time, top):
        """
        Computes for a given time the sea level according to input file parameters.

        Parameters
        ----------
        float : time
            Requested time for which to compute sea level elevation.

        float : top
            Elevation of the core.
        """

        oldsea = self.sealevel
        if self.seafile == None:
            self.sealevel = self.sea0
        else:
            if time < self.seatime.min():
                time = self.seatime.min()
            if time > self.seatime.max():
                time = self.seatime.max()
            self.sealevel = self.seaFunc(time)
        if oldsea == None:
            depth = top
        else:
            depth = top+(self.sealevel-oldsea)

        factors = numpy.ones(self.speciesNb,dtype=float)

        for s in range(self.speciesNb):
            if depth<self.xd[0] and self.edepth[s,1] == self.edepth[s,0]:
                factors[s] = 1.
            elif depth<self.xd[0] and self.edepth[s,1] != self.edepth[s,0]:
                factors[s] = 0.
            elif depth>self.xd[-1] and self.edepth[s,2] == self.edepth[s,3]:
                factors[s] = 1.
            elif depth>self.xd[-1] and self.edepth[s,2] != self.edepth[s,3]:
                factors[s] = 0.
            else:
                factors[s] = self._extract_enviParam( self.xd, self.dtrap[s], depth )

        return depth,factors

    def getSed(self, time):
        """
        Computes for a given time the sediment input according to input file parameters.

        Parameters
        ----------
        float : time
            Requested time for which to compute sediment input.
        """

        if self.sedfile == None:
            self.sedlevel = self.sed0
        else:
            if time < self.sedtime.min():
                time = self.sedtime.min()
            if time > self.sedtime.max():
                time = self.sedtime.max()
            self.sedlevel = self.sedFunc(time)

        factors = numpy.ones(self.speciesNb,dtype=float)
        for s in range(self.speciesNb):
            if self.sedlevel<self.xs[0] and self.esed[s,1] == self.esed[s,0]:
                factors[s] = 1.
            elif self.sedlevel<self.xs[0] and self.esed[s,1] != self.esed[s,0]:
                factors[s] = 0.
            elif self.sedlevel>self.xs[-1] and self.esed[s,2] == self.esed[s,3]:
                factors[s] = 1.
            elif self.sedlevel>self.xs[-1] and self.esed[s,2] != self.esed[s,3]:
                factors[s] = 0.
            else:
                factors[s] = self._extract_enviParam( self.xs, self.strap[s], self.sedlevel )

        return self.sedlevel,factors

    def getFlow(self, time):
        """
        Computes for a given time the flow velocity according to input file parameters.

        Parameters
        ----------
        float : time
            Requested time for which to compute flow velocity value.
        """

        if self.flowfile == None:
            self.flowlevel = self.flow0
        else:
            if time < self.flowtime.min():
                time = self.flowtime.min()
            if time > self.flowtime.max():
                time = self.flowtime.max()
            self.flowlevel = self.flowFunc(time)

        factors = numpy.ones(self.speciesNb,dtype=float)
        for s in range(self.speciesNb):
            if self.flowlevel<self.xf[0] and self.eflow[s,1] == self.eflow[s,0]:
                factors[s] = 1.
            elif self.flowlevel<self.xf[0] and self.eflow[s,1] != self.eflow[s,0]:
                factors[s] = 0.
            elif self.flowlevel>self.xf[-1] and self.eflow[s,2] == self.eflow[s,3]:
                factors[s] = 1.
            elif self.flowlevel>self.xf[-1] and self.eflow[s,2] != self.eflow[s,3]:
                factors[s] = 0.
            else:
                factors[s] = self._extract_enviParam( self.xf, self.ftrap[s], self.flowlevel )

        return factors
