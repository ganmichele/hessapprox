#!/bin/bash

# Tests for Hessian update with Bofill method

echo "******Hessian update******"
echo 
../bin/Bofill_fill_in_H.py  -X cart_water/h2o_cart.xyz  -G cart_water/h2o_grad-cart.dat  -H cart_water/h2o_hess-cart.dat -O cart_water/h2o_approx.dat  -N 20 --xyz

echo "Hessian update done! updated hessians in 'cart_water/h2o_approx.dat'"
