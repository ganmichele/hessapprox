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

- Python3.7+
- numpy (recommended 1.17.4+)
- scipy (optional, recommended 1.3.3+)
- f2py and a reasonably recent gfortran compiler (optional)
- matplotlib (optional, for the optional animation, recommended 3.1.2+)
- ffmpeg (optional, only for the animation, recommended 4.2.4+)

### (Optional) Fortran distance matrix calculator
Most of the code is written in Python, however, optionally, you can
use fortran sources for faster distance calculation. This step is not
mandatory, but could be useful for further developments.

fortran sources for "faster distance calculation" are in ./src.
```bash
cd src
make -f Makefile_calc_dist
make -f Makefile_calc_dist install
```

### (Optional) Distance matrix with scipy
If you do not wish to use the custom distance calculator you should
install scipy with your favourite manager for Python modules.
The scipy and fortran sources are comparable in computational speed.
If none is available the distances would be computed naively and 
inefficiently by a Python function.

### Make the programs available for Python
In ./hessapprox you can find the NGas, DBH and hessian_updt libraries.
The high level programs are in ./bin
You may want to add the hessapprox directory to your PYTHONPATH. For
instance if your shell is bash:
```bash
echo '$PYTHONPATH/the/path/to/hessapprox/' >> ~/.bashrc
source ~/.bashrc
```
or install the modules in a location available to your PYTHONPATH.

## Testing
To test our programs you do not have to add hessapprox to your PYTHONPATH:
```bash
cd tests/
./runme.sh
./runme2.sh
./movie.py
```
runme.sh writes the trajectory files
- cartesian_water/neurons_traj-cart.xyz
- cartesian_water/dbq_traj-cart.xyz
- cartesian_water/h2o_traj.xyz

You may want to visualize them together (for instance with VMD).

runme2.sh writes the approximate hessian matrix in
cart_water/h2o_approx.dat

to run movie.py you need the optional packages scipy and matplotlib.
It generates a movie that shows a NGas optimization for H2O in 2 dimensions.

## Programs usage
In bin you can find the high level programs:
- locate_neurons.py
- locate_DBq.py
- fill_in_H.py

locate_neurons and locate_DBq take as input a space-separated values
file of configurations and output the optimal configurations obtained by
the NGas and DBH methods respectively.
Optionally, if you your configurations are in Cartesian coordinates, you
you can pass the --xyz flag.
Run 
```bash
./locate_neurons.py --help
```
or 
```bash
./locate_DBq.py --help
```
for a thorough list of optional arguments.

fill_in_H takes as input three space-separated values files, containing
respectively: the configurations, gradients and partially filled hessian.
It outputs the hessian matrices as a space separated values file. Run
```bash
./fill_in_H.py --help
```
for a thorough list of the arguments
