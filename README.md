# fabsfit
A small program to fit a parametrization of the the absorptive scattering factor.


## Credits

The workflow of this code was designed by [Dr. José Ángel Castellanos-Reyes](https://www.uu.se/en/contact-and-organisation/staff?query=N22-1026). We acknowledge the [abTEM](https://github.com/abTEM/abTEM) multislice code for the compilation of parametrization data of elastic scattering factors.


## Installation

The installation can be performed using pip:
```
$ pip install fabsfit
```


### Elastic Parametrization Data Download

We use the parametrization data of [abTEM](https://github.com/abTEM/abTEM/tree/main/abtem/parametrizations/data). To download the files to the fabsfit data directory, run
```
$ fabsfit download
```
Alternatively you can save the data files also in any other directory using
```
$ fabsfit download /path/to/dir
```


## Usage

fabsfit can be run in two ways: directly from the command line or as a Python package. Below we give examples of the syntax for both cases.


### Command line interface

We can perform a fit of the absorptive scattering factor for Manganese (Mn) and a isotropic $B_{\mathrm{iso}}$ parameter of the Debye-Waller factor of 0.39 Å<sup>2 </sup> at the default electron energy of 100 keV with the command
```
$ fabsfit fit Mn 0.39
```
Help for more advanced usage options can be obtained with
```
$ fabsfit fit --help
```


### Python

The functionality in this package can be directly exposed in Python using
```
from fabsfit.parametrization import Parametrization

# Initializae a Parametrization object
param = Parametrization('doyleturner', 'doyleturner', 'Mn', 100, 0.39)

# perform the fit
param.fit()
```
The resulting fit would be equivalent to the fit in the example for the command line interface above.


## Background

The total scattering factor can be written as
```math
f_{\mathrm{tot}} (s) = \left[ f_{\mathrm{el}}(s) + i f_{\mathrm{abs}} (s) \right] \exp(-2 B_{\mathrm{iso}} s^2)
```
for the purpose of electron scattering simulations, where $f_{\mathrm{el}}$ is the elastic scattering factor and $f_{\mathrm{abs}}$ is the absorptive scattering factor. Different inealstic scattering processes contribute to the absorptive scattering factor, such as core-electron excitations, plasmon excitations, and phonon (vibrational) excitations.

Here we focus on vibrational excitations and calculate the absorptive scattering factor for uncorrelated atomic motion as first considered by Hall and Hirsch (Hall and Hirsch, [1965](https://doi.org/10.1098/rspa.1965.0136)). 
```math
f_{\mathrm{abs}} (s) = \frac{4 \pi \hbar}{m_0 V} \int_{\Omega} d^2 s' \, f_{\mathrm{el}} \left( \left| \frac{s}{2} + s' \right| \right) 
f_{\mathrm{el}} \left( \left| \frac{s}{2} - s' \right| \right) \left\{ 1 - \exp \left[ -2B_{\mathrm{iso}} \left( s'^2 - \frac{s^2}{4} \right) \right] \right\}
```


Right now only the *Doyle-Turner* form, i.e.,
```math
f (s) = \sum_{i=1}^{5} a_i \exp(- b_i s^2)
```
is supported, both for the elastic scattering factor $f_{\mathrm{el}}$ and the absorptive scattering factor $f_{\mathrm{abs}}$. We plan to extend this in the future.





## Development

For development purposes, clone the repository
```
git clone https://github.com/pzeiger/fabsfit.git
```
and build the python package wheel using [hatch](https://hatch.pypa.io)
```
hatch -c build
```






