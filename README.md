# Hessian approximation methods for molecular dynamics simulations
this code implements:</br>
the NGas method introduced in ...</br>
the DBH method introduced in ...</br>
the hessian update methods ...</br>

----

# Table of Contents
1. [Installation requirements](#Installation-requirements)
2. [Testing](#Testing)
3. [Programs usage](#Programs-usage)

----

## Installation requirements
To run the codes you need:

- Python3.6+
- numpy 
- (scipy)
- (f2py and a reasonably recent gfortran compiler)

### (Optional) Fortran distance matrix calculator
Most of the code is written in Python, however, optionally, you can
use fortran sources for faster distance calculation. This step is not
mandatory, but could be useful for further developments.

fortran sources for "faster distance calculation" are in ./src.
```bash
cd src
make -f Makefile_calc_dist
make install
```

### (Optional) Distance matrix with scipy
If you do not wish to use the custom distance calculator you should
install scipy with your favourite manager for Python modules.
The scipy and fortran sources are comparable in computation speed.
If none is available the distances would be computed naively and 
inefficiently by a Python function.

### Make the programs available for Python
In ./hessapprox you can find the NGas, DBH and hessian_updt libraries.
The high level programs are in ./bin
You may want to add the hessapprox directory to your PYTHONPATH. For
instance if your shell is bash:
```bash
echo '$PYTHONPATH/the/path/to/hessapprox/hessapprox/' >> ~/.bashrc
source ~/.bashrc
```
or install the modules in a location available to your PYTHONPATH.

## Testing
To test our programs you do not have to add hessapprox to your PYTHONPATH:
```bash
cd tests/
./runme.sh
./runme2.sh
```
runme.sh writes the trajectory files
- cartesian_water/neurons_traj-cart.xyz
- cartesian_water/dbq_traj-cart.xyz
- cartesian_water/h2o_traj.xyz

You may want to visualize them together (for instance with VMD).

## Programs usage
Further tests can be done straightforwardly:

locate_neurons and locate_DBq take as input a space-separated values
file of configurations and output the optimal configurations obtained by
the NGas and DBH methods respectively.
For a thorough list of optional arguments run
```bash
./locate_neurons.py --help
```
or 
```bash
./locate_DBq.py --help
```

fill_in_H takes as input three space-separated values files, containing
respectively: the configurations, gradients and partially filled hessian.
It outputs the hessian matrices as a space separated values file. Run
```bash
./fill_in_H.py -X q.dat -G g.dat -H H.dat [-N num evals] [-L lambda] [--xyz]
```
