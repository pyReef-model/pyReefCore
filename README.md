# pyReefCore
Synthetic coral reef core model

# Generalized Lotka-Volterra model

The most common models of species evolution in ecological modeling are the predator-prey **Lotka-Volterra (LV)** equation and its modifications.

From LV equations, one can formulate the **Generalized Lotka-Volterra (GLV)** equation that allows unlimited number of species and different types of interactions among species. The GLV equation is mainly formed by two parts, the logistic growth/decay of a species and its interaction with the other species,

![equation](http://bit.ly/2mgyrpI)
$$\frac{dx_i}{dt} = \epsilon_i x_i + \sum_{j=1}^{N_s} \alpha_{ij}x_ix_j$$

where $x_i$ is the population density of species _i_; $\epsilon_i$ is the intrinsic rate of increase/decrease of a population of species _i_ (also called **Malthusian** parameter); $\alpha_{ij}$ is the interaction coefficient among the species association _i_ and _j_, (a particular case is $\alpha_{ii}$, the interaction of one species association with itself); and _t_ is time. Equation 1 can be written in matrix formulation as:

$$\frac{dx_i}{dt} = diag[X](\epsilon + AX)$$

where $X$ is the vector of population densities of each species _i_; $\epsilon$ is the vector of all _Mathusian_ parameters; $A$ is the matrix of interaction coefficients, also known as community matrix; and $diag[X]$ is a square matrix with diagonal elements equal to $X$, and zeros outside the diagonal.

Here we show how to solve this system of ODEs with 5 species.

## Definition of species rate and community matrix

To solve the ODEs we first need to define some initial conditions:

- the initial species population number $X0$
- the intrinsic rate of a population species $\epsilon$
- the interaction coefficients among the species association $\alpha$

```python
%matplotlib inline
import matplotlib.pyplot as plt

import numpy
#import cmocean

%config InlineBackend.figure_format = 'svg'
from pyReefCore.model import Model

# initialise model
reef = Model()

# Define the XmL input file
reef.load_xml('input.xml')

# Run to a given time
reef.run_to_time(500,showtime=100.)
```

## Solving the ODEs system

The mathematical model for the species population evolution results in a set of differential equations (ODEs), one for each species associations modeled. The **Runge-Kutta-Fehlberg** method (_RKF45_ or _Fehlberg_ as defined in the odespy library) is used to solve this **GLV ODE system**.

The _Fehlberg_ method requires 5 parameters to solve the GLV ODEs:

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

We define the maximum carbonate production rate (m/y) for each species
