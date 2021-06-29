#!/bin/bash

# Collection of tests

# Normal modes coordinates
echo "******NGas FROM NORMAL MODES******"
echo 
../bin/locate_neurons.py nm_water/h2o_NMq.dat -N 100 -O nm_water/neurons -T -i 0:1000
echo
echo "Converting neurons to Cartesian. Output in cart_water/neurons_traj-cart.xyz"
echo 

../bin/toCart.py nm_water/neurons-traj.dat  nm_water/cnorm.dat  nm_water/mass.dat  nm_water/qrt.dat  cart_water/neurons_traj-cart.xyz

echo "******DBH FROM NORMAL MODES******"
echo 
../bin/locate_DBq.py nm_water/h2o_NMq.dat -R 2.98 -O nm_water/dbq -T -i 0:1000
echo
echo "Converting neurons to Cartesian. Output in cart_water/dbq_traj-cart.xyz"
echo 

../bin/toCart.py  nm_water/dbq  nm_water/cnorm.dat  nm_water/mass.dat  nm_water/qrt.dat  cart_water/dbq_traj-cart.xyz

# NGas with Cartesian coordinates (not recommanded but works ok)
echo
echo "******NGas FROM CARTESIAN (RECOMMANDED FOR TESTING PURPOSES ONLY)******"
echo
../bin/locate_neurons.py cart_water/h2o_traj.xyz --xyz -S cent -N 100 -O cart_water/neurons-cart -T

echo 
echo "You may compare cart_water/neurons_traj-cart.xyz, cart_water/dbq_traj-cart.xyz and cart_water/h2o_traj.xyz"
echo 

