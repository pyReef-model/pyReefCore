##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the pyReefCore synthetic coral reef core model app.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
   pyReefCore Model main entry file.
"""
import time
import numpy as np
import mpi4py.MPI as mpi

from pyReefCore import (preProc, xmlParser, enviForce, coralGLV, coreData, modelPlot)

# profiling support
import cProfile
import os
import pstats
import StringIO


class Model(object):
    """State object for the pyReef model."""

    def __init__(self):
        """
        Constructor.
        """

        # Simulation state
        self.dt = 0.
        self.tNow = 0.
        self.tDisp = 0.
        self.waveID = 0
        self.outputStep = 0
        self.applyDisp = False
        self.simStarted = False

        self.dispRate = None

        self._rank = mpi.COMM_WORLD.rank
        self._size = mpi.COMM_WORLD.size
        self._comm = mpi.COMM_WORLD

        # Initialise pre-processing functions
        self.enviforcing = preProc.preProc()

    def load_xml(self, filename, verbose=False):
        """
        Load an XML configuration file.
        """

        # Only the first node should create a unique output dir
        self.input = xmlParser.xmlParser(filename, makeUniqueOutputDir=(self._rank == 0))
        self.tNow = self.input.tStart
        self.tCoral = self.tNow
        self.tLayer = self.tNow #+ self.input.laytime

        # Seed the random number generator consistently on all nodes
        seed = None
        if self._rank == 0:
            # limit to max uint32
            seed = np.random.mtrand.RandomState().tomaxint() % 0xFFFFFFFF
        seed = self._comm.bcast(seed, root=0)
        np.random.seed(seed)
        self.iter = 0
        self.layID = 0

        # Initialise environmental forcing conditions
        self.force = enviForce.enviForce(input=self.input)

        # Initialise core data
        self.core = coreData.coreData(input=self.input)
        # Environmental forces functions
        self.core.seatime = self.force.seatime
        self.core.sedtime = self.force.sedtime
        self.core.flowtime = self.force.flowtime
        self.core.seaFunc = self.force.seaFunc
        self.core.sedFunc = self.force.sedFunc
        self.core.flowFunc = self.force.flowFunc

        # Initialise plotting functions
        self.plot = modelPlot.modelPlot(input=self.input)

        return

    def run_to_time(self, tEnd, showtime=10, profile=False, verbose=False):
        """
        Run the simulation to a specified point in time (tEnd).

        If profile is True, dump cProfile output to /tmp.
        """

        timeVerbose = self.tNow+showtime

        if profile:
            pid = os.getpid()
            pr = cProfile.Profile()
            pr.enable()

        if self._rank == 0:
            print 'tNow = %s [yr]' %self.tNow

        if tEnd > self.input.tEnd:
            tEnd = self.input.tEnd
            print 'Requested time is set to the simulation end time as defined in the XmL input file'

        if self.tNow == self.input.tStart:
            # Initialise Generalized Lotka-Volterra equation
            self.coral = coralGLV.coralGLV(input=self.input)
            self.odeRKF = self.coral.solverGLV()

        # Perform main simulation loop
        # NOTE: number of iteration for the ODE during a given time step, could be user defined...
        N = 100

        # Define environmental factors
        dfac = np.ones(self.input.speciesNb,dtype=float)
        sfac = np.ones(self.input.speciesNb,dtype=float)
        ffac = np.ones(self.input.speciesNb,dtype=float)
        while self.tNow < tEnd:

            # Initial coral population
            if self.tNow == self.input.tStart:
                self.coral.population[:,self.iter] = self.input.speciesPopulation

            # Get sea-level
            if self.input.seaOn:
                tmp = self.core.topH
                self.core.topH, dfac = self.force.getSea(self.tNow, tmp)

            # Get sediment input
            if self.input.sedOn:
                sfac = self.force.getSed(self.tNow)

            # Get flow velocity
            if self.input.flowOn:
                ffac = self.force.getFlow(self.tNow)

            # Limit species activity from environmental forces
            tmp = np.minimum(dfac, sfac)
            fac = np.minimum(ffac, tmp)
            self.coral.epsilon = self.input.malthusParam*fac

            # Initialise RKF conditions
            self.odeRKF.set_initial_condition(self.coral.population[:,self.iter])

            # Define coral evolution time interval and time stepping
            self.tCoral += self.input.tCarb
            tODE = np.linspace(self.tNow, self.tCoral, N+1)
            self.dt = tODE[1]-tODE[0]

            # Solve the Generalized Lotka-Volterra equation
            coral,t = self.odeRKF.solve(tODE)
            population = coral.T

            # Update coral population
            self.iter += 1
            ids = np.where(self.coral.epsilon==0.)[0]
            population[ids,-1] = 0.
            ids = np.where(np.logical_and(fac==1,population[:,-1]==0.))[0]
            population[ids,-1] = 1.
            self.coral.population[:,self.iter] = population[:,-1]

            # Compute carbonate production and update coral core characteristics
            self.core.coralProduction(self.layID,self.coral.population[:,self.iter],
                                      self.coral.epsilon)

            # Update time step
            self.tNow = self.tCoral

            # Update stratigraphic layer ID
            if self.tLayer < self.tNow :
                self.tLayer += self.input.laytime
                self.layID += 1

            if self._rank == 0 and self.tNow>=timeVerbose:
                timeVerbose = self.tNow+showtime
                print 'tNow = %s [yr]' %self.tNow

        # Update plotting parameters
        self.plot.pop = self.coral.population
        self.plot.timeCarb = self.coral.iterationTime
        self.plot.depth = self.core.thickness
        self.plot.sedH = self.core.coralH
        self.plot.timeLay = self.core.layTime
        self.plot.surf = self.core.topH

        return

    def ncpus(self):
        """
        Return the number of CPUs used to generate the results.
        """

        return 1
