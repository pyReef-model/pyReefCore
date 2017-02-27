##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the pyReefCore synthetic coral reef core model app.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Functions used to create environmental forcing conditions for pyReefCore simulation.
"""

import errno
import matplotlib
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

class preProc:
    """
    Class for creating pyReefCore environmental forcing conditions.
    """

    def __init__(self, curve=None):
        """
        Initialization function which could take an existing curve as an option. If you build
        your own sea-level curve you need to define this parameter as None.

        Parameters
        ----------
        variable : curve
            An existing curve for sea-level or flow velocity or sedimentation rate.
            The file is defined with 2 columns:
                + Column 1 = Time (a)
                + Column 2 = Curve value for the considered time [m or m/d]
        """

        self.df = None
        self.time = None
        self.func = None

        if curve != None:
            self.build = False
            self.df = pd.read_csv(curve1, sep=r'\s+', header=None, names=['h','t'])
        else:
            self.build = True
            self.func = None
            self.time = None

        return

    def buildCurve(self, timeExt = None, timeStep = None, funcExt = None,
                   ampExt = None, periodExt = None):
        """
        Curve created which interpolate linearly the averaged values of the environmental
        parameter trends over the specified time period.

        Parameters
        ----------
        variable: timeExt
            Extent of the simulation time: start/end time (in years)

        variable: timeStep
            Discretisation step for time range (in years).

        variable: funcExt
            Environmental factor value for starting and ending times (in metres)

        variable: ampExt
            Amplitude of the environmental factor wave for starting and ending times (in metres)

        variable: periodExt
            Period of the nvironmental factor wave for starting and ending times (in years)
        """

        dt = float(timeStep)
        so = float(funcExt[0])
        sm = float(funcExt[1])
        to = float(timeExt[0])
        tm = float(timeExt[1])
        Ao = float(ampExt[0])
        Am = float(ampExt[1])
        Po = float(periodExt[0])
        Pm = float(periodExt[1])

        self.time = np.arange(to,tm+dt,dt,dtype=np.float)

        # Environmental factor
        a0 = (sm - so)/(tm - to)
        b0 = so - a0 * to
        self.func = a0 * self.time + b0
        # Amplitude
        a1 = (Am - Ao)/(tm - to)
        b1 = Ao - a1 * to
        A = a1 * self.time + b1
        # Period
        a2 = (Pm - Po)/(tm - to)
        b2 = Po - a2 * to
        P = a2 * self.time + b2
        # Enveloppe
        self.env1 = a0 * self.time + b0 - 1.
        self.env2 = a0 * self.time + b0 + 1.

        for t in range(len(self.time)):
            self.func[t] += A[t] * np.cos(2.* np.pi * (self.time[t] - to) / P[t])


        f = interpolate.interp1d(self.time, self.func, kind='cubic')
        tnew = np.arange(self.time.min(),timeExt[1]+dt/10.,dt/10.,dtype=np.float)
        self.time = tnew
        self.func = f(tnew)

        return

    def readCurve(self, timeStart = None, timeEnd = None, dt = 10.):
        """
        Read environmental function curves.

        Parameters
        ----------
        variable: timeStart
            Simulation time start in years.

        variable: timeEnd
            Simulation time end in years.

        variable: dt
            Discretisation step for time range (in Ma).
        """

        self.df.columns[1]
        list1 = list(self.df)
        time = self.df[list1[0:len(list1)-1]].values[:,0]
        func = self.df[list1[len(list1)-1]].values

        if timeStart == None:
            timeStart = time.min()
        if timeEnd == None:
            timeEnd = time.max()

        interpFn = interpolate.interp1d(self.time, self.func)

        self.time = np.arange(timeStart, timeEnd+dt, dt)
        self.func = interpFn(self.time)

        return

    def plotCurves(self, size=(5,5), lwidth = 3, title=None, color=None, font=9, dpi=80, figName = None):
        """
        Plot environmental curves.

        Parameters
        ----------
        variable: size
            Size of the figure to plot.

        variable : lwidth
            Figure lines width

        variable : title
            Figure title

        variable: color
            Matplotlib color to use.

        variable : font
            Figure font size

        variable : dpi
            Figure resolution

        variable: figName
            Name of the saved file.
        """

        matplotlib.rcParams.update({'font.size': font})

        # Define figure size
        fig, ax = plt.subplots(1,figsize=size, dpi=dpi)
        ax.set_axis_bgcolor('#f2f2f3')
        fig.tight_layout()

        # Plotting curve
        plt.plot(self.func, self.time, color=color, linewidth=lwidth)

        if title != None:
            titlepos = plt.title(title, fontsize=font+3,fontweight='bold')
            titlepos.set_y(1.02)

        plt.xlabel('Environmental factor',fontsize=font+3)
        plt.ylabel('Time [a]',fontsize=font+3)
        ax.set_ylim(self.time[0], self.time[-2])
        plt.grid()
        plt.tick_params(axis='both', which='major', labelsize=font)
        plt.show()

        if figName is not None:
            fig.savefig(figName, dpi = dpi)

        return

    def exportCurve(self, factor=1., nameCSV='function.csv'):
        """
        Write CSV sea level file following pyReef core requirements:
            + 2 columns file containing time in years (1st column) and environmental parameter (2nd column),
            + time will be in increasing order starting at the oldest time,
            + past times are negative,
            + the separator is a space.

        Parameters
        ----------
        variable: curve
            Environmental parameter to save.

        variable : factor
            Factor to convert from given time unit to years (ex: Ma -> a).

        variable: nameCSV
            Name of the saved CSV file.
        """

        df = pd.DataFrame({'X':np.around(self.time*factor, decimals=0),'Y':np.around(self.func, decimals=3)})
        df.to_csv(str(nameCSV),columns=['X', 'Y'], sep=' ', index=False ,header=0)

        return
