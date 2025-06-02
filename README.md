# **RECOVER** 
**RE**treving **C**hemical **O**rdering **V**ia **E**xplicating **R**adial-distribution function.  <br> 
Package for retreving ***Warren-Cowley SRO parameters*** from perturbed ***atom-probe tomography (APT) data.*** The main funcitons of this package are:

* Compute fkr functions through simulating perturbed APT data through either isotropic or anisotropic perturbation of atomic positions. (_examples/numerical-simulation_)
* Retreive SRO parameters from experimental APT data using the calculated fkr functions (_examples/experimental-APT-analysis_)

Theoretical foundation and results of this package is documented at (URL of this paper). 
  
## Requirements:
  * Python >= 3.10
  * Cython
  * Numpy
  * Scipy >=1.15.2


## Installation:
  From inside of the RECOVER root directory, run the following python command.
  * python setup.py install

## Usage:
  * Examples documenting the package usage are located in the examples directory.

## Examples:
  * numerical-simulation: simulate raw APT data and compute fkr functions through perturbation of GeSn configurations
  * experimental-APT-analysis: retrieve SRO parameters from experimental APT raw data of GeSn alloy based on the calculated fkr functions
