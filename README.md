# Hessian approximation methods for molecular dynamics simulations
this code implements:</br>
- the NGas method introduced in Gandolfi M, Ceotto M J. Chem. Theory Comput. 153, 204104 (2021)</br>
- the DBH method introduced in Conte R, Gabas F, Botti G, Zhuang Y, Ceotto M J. Chem. Phys. 150, 244118 (2019)</br>
- the hessian update methods described in H Wu, M Rahman, J Wang, U Louderaj, W L Hase, and Y Zhuang. J. Chem. Phys. 133, 074101 (2010) 

**If you use the codes provided here or parts of them, please cite the following
article:**</br>
Gandolfi M, Ceotto M J. Chem. Theory Comput. (2021)</br>

----

# Table of Contents
1. [Installation requirements](#Installation-requirements)
2. [Testing](#Testing)
3. [Programs usage](#Programs-usage)
4. [Authors](#Authors)

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
If you do not wish to use our Fortran distance calculator you should
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
echo "export PYTHONPATH=${PYTHONPATH}:/the/path/to/hessapprox/" >> ~/.bashrc
source ~/.bashrc
```
or install the modules in a location available to your PYTHONPATH.

## Testing
To test our programs you do not have to add hessapprox to your PYTHONPATH:
```bash
cd tests/
./runme.sh
./runme2.sh
./runme3.sh
./movie.py
```
runme.sh writes the trajectory files
- cartesian_water/neurons_traj-cart.xyz
- cartesian_water/dbq_traj-cart.xyz
- cartesian_water/h2o_traj.xyz

You may want to visualize them together (for instance with VMD).

runme2.sh writes the approximate hessian matrix in
cart_water/h2o_approx.dat

runme3.sh tests the NGas method on Cartesian inputs.

To run movie.py you need the optional packages scipy and matplotlib.
It generates a short movie that shows a NGas optimization for H2O in 2 
dimensions.

## Programs usage
In bin you can find the high level programs:
- locate_NGas.py
- locate_DBq.py
- Bofill_fill_in_H.py

locate_NGas and locate_DBq take as input a space-separated values
file of configurations. They output the optimal configurations obtained by
the NGas and DBH methods respectively, as well as a file (called
relations_NGas.dat / relations_DBq.dat) that says which trajectory geometry
belongs to which neuron / DBq.
Optionally, if you your input configurations are in Cartesian coordinates,
you you can pass the --xyz flag.
Run 
```bash
./locate_NGas.py --help
```
or 
```bash
./locate_DBq.py --help
```
for a thorough list of optional arguments.

Bofill_fill_in_H takes as input three space-separated values files, containing
respectively: the configurations, gradients and incomplete list of hessians.
It outputs the hessian matrices as a space separated values file. Run
```bash
./Bofill_fill_in_H.py --help
```
for a thorough list of the arguments

## Authors

* **Michele Gandolfi** 
* **Michele Ceotto**

For more information about the authors and the research group visit our website
[Ceotto Group](https://sites.unimi.it/ceotto/)

