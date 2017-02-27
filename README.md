# pyReefCore

<div align="center">
    <img width=600 src="http://sydney.edu.au/science/geosciences/images/core.jpg" alt="OTR core from the Geocoastal Group - USyD" title="Core sample, including a well-preserved Faviidae coral recovered from 16 m depth. The red arrows are drawn on to indicate upwards recovery direction."</img>
</div>


**pyReefCore** is a 1D model which simulates evolution of mixed carbonate-siliciclastic system under environmental forcing conditions (_e.g._ sea-level, water flow, siliciclastic input). The carbonate production model simulates the logistic growth and interaction among species based on **Generalized Lotka-Volterra** equations. The environmental forces are converted to factors and combined into one single environmental value to model the evolution of species. The interaction among species is quantified using a _community matrix_ that captures the beneficial or detrimental effects of the presence of each species on the other.

## Generalized Lotka-Volterra model

The most common models of species evolution in ecological modeling are the predator-prey **Lotka-Volterra (LV)** equation and its modifications.

From LV equations, one can formulate the **Generalized Lotka-Volterra (GLV)** equation that allows unlimited number of species and different types of interactions among species. The GLV equation is mainly formed by two parts, the logistic growth/decay of a species and its interaction with the other species,

<p align="center"><img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/36c0f67b25b915d50ae129824705fc29.svg?invert_in_darkmode" align=middle width=173.0751pt height=50.188545pt/></p>

where <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/9fc20fb1d3825674c6a279cb0d5ca636.svg?invert_in_darkmode" align=middle width=13.993485pt height=14.10255pt/> is the population density of species _i_; <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/1cd32b0756da515bc59142b9318ff797.svg?invert_in_darkmode" align=middle width=11.28105pt height=14.10255pt/> is the intrinsic rate of increase/decrease of a population of species _i_ (also called **Malthusian** parameter); <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/8175b4b012861c57d7f99a503fdcaa72.svg?invert_in_darkmode" align=middle width=21.19161pt height=14.10255pt/> is the interaction coefficient among the species association _i_ and _j_, (a particular case is <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/ff41937f5cd113c5b9d670fd51ac28f1.svg?invert_in_darkmode" align=middle width=19.743405pt height=14.10255pt/>, the interaction of one species association with itself); and _t_ is time. Equation 1 can be written in matrix formulation as:

<p align="center"><img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/75aeba56acecbc8c3ead4fd18a21d6b3.svg?invert_in_darkmode" align=middle width=169.01445pt height=33.769395pt/></p>

where <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.85297pt height=22.38192pt/> is the vector of population densities of each species _i_, <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.6475035pt height=14.10255pt/> is the vector of all _Mathusian_ parameters, <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode" align=middle width=12.282765pt height=22.38192pt/> is the matrix of interaction coefficients, also known as community matrix, and <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/9f0dfe8a0e850780e96018103aa5fe64.svg?invert_in_darkmode" align=middle width=55.17996pt height=24.56553pt/> is a square matrix with diagonal elements equal to <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/cbfb1b2a33b28eab8a3e59464768e810.svg?invert_in_darkmode" align=middle width=14.85297pt height=22.38192pt/>, and zeros outside the diagonal.

## Definition of species rate and community matrix

To solve the ODEs, the user needs to define several initial conditions:

- the initial species population number <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/8058174c3e47972feecfee6a81720995.svg?invert_in_darkmode" align=middle width=23.046375pt height=22.38192pt/>
- the intrinsic rate of a population species <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/7ccca27b5ccc533a2dd72dc6fa28ed84.svg?invert_in_darkmode" align=middle width=6.6475035pt height=14.10255pt/>
- the interaction coefficients among the species association <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode" align=middle width=10.537065pt height=14.10255pt/>

Several other input are required and will need to be set in the **XmL** inputfile. An example of such file is provided [here](https://github.com/pyReef-model/pyReefCore/blob/master/Tests/input.xml).

## Solving the ODEs system

The mathematical model for the species population evolution results in a set of differential equations (ODEs), one for each species associations modeled. The **Runge-Kutta-Fehlberg** method (_RKF45_ or _Fehlberg_ as defined in the [**odespy**](http://hplgit.github.io/odespy/doc/pub/tutorial/html/main_odespy.html) library) is used to solve the **GLV ODE system**.

The _Fehlberg_ method requires 5 parameters :

- the step-size (or time step)
- an initial population of the species association
- the relative tolerance for solution
- the absolute tolerance for solution
- the minimum step size for an adaptive algorithm.

## Carbonate production

Once a species association population is computed, carbonate production is calculated using a carbonate production factor. Production factors are specified for the maximum population, and linearly scaled to the actual population following the relation
<p align="center"><img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/86747a51e4fe91aa94ce20bfe5b4a600.svg?invert_in_darkmode" align=middle width=106.14813pt height=36.235155pt/></p>

where <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/df5a289587a2f0247a5b97c1e8ac58ca.svg?invert_in_darkmode" align=middle width=12.78882pt height=22.38192pt/> is the carbonate production, <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.913963pt height=20.1465pt/> is time, <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/ea449f9e9a48e2959872aac8fa65e1ca.svg?invert_in_darkmode" align=middle width=38.586405pt height=22.38192pt/> is the carbonate production factor when population is at its maximum, and <img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/655ca15e2b101fb431577b12d4442580.svg?invert_in_darkmode" align=middle width=18.5427pt height=22.38192pt/> is the maximum population of species _i_, computed as

<p align="center"><img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/43e99ce899b3ecaaa637f00ee554b2ae.svg?invert_in_darkmode" align=middle width=63.86358pt height=31.913475pt/></p>

which gives:

<p align="center"><img src="https://rawgit.com/pyReef-model/pyReefCore/master/svgs/2baed347518e8396d33d82d9c18e283b.svg?invert_in_darkmode" align=middle width=124.960935pt height=36.235155pt/></p>

We define the maximum carbonate production rate (m/y) for each species in the **XmL** input file.

## Installation

### Local install

The code is available from our github [page](https://github.com/pyReef-model/pyReefCore.git) and can be obtained either from this page or using **git**
```
git clone https://github.com/pyReef-model/pyReefCore.git
```

Once donwloaded, navigate to the **pyReefCore** folder and run the following command:
```
pip install -e /workspace/volume/pyReefCore/
```

### Docker container

The code is available from Docker Hub at [pyreefmodel/pyreef-Docker](https://hub.docker.com/u/pyreefmodel/) and can be downloaded using **Kitematic**. An example of model is provided in the [Tests](https://github.com/pyReef-model/pyReefCore/tree/master/Tests) folder using IPython Notebook.

## Usage

pyReefCore can be use from an _IPython notebook_ or a _python script_ directly. An example of functions available is provided below:

```python

%matplotlib inline
import matplotlib.pyplot as plt

import numpy

%config InlineBackend.figure_format = 'svg'
from pyReefCore.model import Model

# Initialise model
reef = Model()

# Define the XmL input file
reef.load_xml('input.xml')

# Run to a given time (for example 500 years)
reef.run_to_time(500.,showtime=100.)

# Define a colorscale to display the core
# Some colormaps are available from the following link:
# http://matplotlib.org/examples/color/colormaps_reference.html
from matplotlib.cm import terrain
nbcolors = len(reef.core.coralH)+3
colors = terrain(numpy.linspace(0, 1, nbcolors))

# Plot evolution of species population with time
reef.plot.speciesTime(pop=reef.coral.population, time=reef.coral.iterationTime, colors=colors,
                      size=(8,4), font=8, dpi=80,fname='pop.pdf')

# Plot evolution of species population with depth
reef.plot.speciesDepth(pop=reef.coral.population, depth=reef.core.thickness, colors=colors,
                       size=(8,4), font=8, dpi=80)

# Plot coral facies distribution core and assemblages
reef.plot.drawCore(pop=reef.core.coralH, depth=reef.core.thickness, depthext = [5.,10.],
                   thext = [0.,0.004], propext = [0.,0.5], surf=reef.core.topH, lwidth = 3,
                   colors=colors, size=(8,10), font=8, dpi=380, figname='out.pdf',
                   filename='out.csv', sep='\t')
```

## Model results
