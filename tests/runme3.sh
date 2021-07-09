#!/bin/bash

# test

# NGas with Cartesian coordinates (not recommanded but works ok)
echo
echo "******NGas FROM CARTESIAN (FOR TESTING PURPOSES ONLY)******"
echo
../bin/locate_NGas.py cart_water/h2o_traj.xyz --xyz -S cent -N 100 -O cart_water/neurons-cart -T

