# pyReefCore

**pyReefCore** is a 1D model which simulates evolution of mixed carbonate-siliciclastic system under environmental forcing conditions (_e.g._ sea-level, water flow, siliciclastic input). The carbonate production model simulates the logistic growth and interaction among species based on **Generalized Lotka-Volterra** equations. The environmental forces are converted to factors and combined into one single environmental value to model the evolution of species. The interaction among species is quantified using a _community matrix_ that captures the beneficial or detrimental effects of the presence of each species on the other.

## Generalized Lotka-Volterra model

The most common models of species evolution in ecological modeling are the predator-prey **Lotka-Volterra (LV)** equation and its modifications.

From LV equations, one can formulate the **Generalized Lotka-Volterra (GLV)** equation that allows unlimited number of species and different types of interactions among species. The GLV equation is mainly formed by two parts, the logistic growth/decay of a species and its interaction with the other species,

$$\frac{dx_i}{dt} = \epsilon_i x_i + \sum_{j=1}^{N_s} \alpha_{ij}x_ix_j$$

where $x_i$ is the population density of species _i_; $\epsilon_i$ is the intrinsic rate of increase/decrease of a population of species _i_ (also called **Malthusian** parameter); $\alpha_{ij}$ is the interaction coefficient among the species association _i_ and _j_, (a particular case is $\alpha_{ii}$, the interaction of one species association with itself); and _t_ is time. Equation 1 can be written in matrix formulation as:

$$\frac{dx_i}{dt} = diag[X](\epsilon + AX)$$

where $X$ is the vector of population densities of each species _i_, $\epsilon$ is the vector of all _Mathusian_ parameters, $A$ is the matrix of interaction coefficients, also known as community matrix, and $diag[X]$ is a square matrix with diagonal elements equal to $X$, and zeros outside the diagonal.

## Definition of species rate and community matrix

To solve the ODEs, the user needs to define several initial conditions:

- the initial species population number $X0$
- the intrinsic rate of a population species $\epsilon$
- the interaction coefficients among the species association $\alpha$

Several other input are required and will need to be set in a **XmL** inputfile. An example of such file is provided in [here](https://github.com/pyReef-model/pyReefCore/blob/master/Tests/input.xml).

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
$$ \frac{dP}{dt} = R_{max}\frac{x_i}{K_i}$$

where $P$ is the carbonate production, $t$ is time, $R_{max}$ is the carbonate production factor when population is at its maximum, and $K_i$ is the maximum population of species _i_, computed as

$$K_i=\frac{\epsilon_i}{\alpha_{ii}}$$

which gives:

$$ \frac{dP}{dt} = R_{max}\frac{\alpha_{ii}\, x_i}{\epsilon_i}$$

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

# initialise model
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
reef.plot.drawCore(pop=reef.core.coralH, depth=reef.core.thickness, surf=reef.core.topH,
                   colors=colors, size=(8,10), font=8, dpi=380, fname='out.pdf')
```
