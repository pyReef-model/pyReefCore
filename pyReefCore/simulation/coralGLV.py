##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the pyReefCore synthetic coral reef core model app.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module solves the Generalized Lotka-Volterra (GLV) equation that allows unlimited
number of species and different types of interactions among species.
The GLV equation is mainly formed by two parts, the logistic growth/decay of a species
and its interaction with the other species.
"""
import os
import numpy
import odespy

class coralGLV:
    """
    This class solves the Generalized Lotka-Volterra equation using Runge-Kutta-Fehlberg
    method (RKF45 or Fehlberg as defined in the odespy library)
    """

    def __init__(self, input = None):
        """
        Constructor.
        """

        # RKF relative tolerance for solution
        self.rtol = 1.e-6
        # RKF absolute tolerance for solution
        self.atol = self.rtol
        # RKF minimum step size for an adaptive algorithm.
        self.min_step = 1.e-4
        # Definition of the intrinsic rate of a population species
        self.epsilon = input.malthusParam
        # Community matrix representing the interactions between species
        self.alpha = input.communityMatrix
        # Coral population record through time
        self.iterationTime = numpy.arange(input.tStart, input.tEnd+input.tCarb, input.tCarb)
        self.population = numpy.zeros((input.speciesNb,len(self.iterationTime)),dtype=float)

        return

    def _functionGLV(self, X, t):
        """
        This function solves the ODEs defining for the Generalized Lotka-Volterra equation.

        Parameters
        ----------

        variable : X
            Species population distribution at current time step.

        variable : t
            Time step on which to solve the ODEs for.
        """

        function = numpy.zeros(len(self.epsilon))

        for eq in range(len(self.epsilon)):
            function[eq] = (self.epsilon[eq]+numpy.sum(self.alpha[eq,:]*X))*X[eq]

        return function

    def solverGLV(self):
        """
        This function build the RKF solver used for the Generalized Lotka-Volterra equation.
        """

        # RKF initialisation
        odeRKF = odespy.Fehlberg(self._functionGLV, atol=self.atol,
                                   rtol=self.rtol, min_step=self.min_step)

        return odeRKF
