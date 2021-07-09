#!/usr/bin/env python
"""
Hessian Approximation Methods. Full notice in LICENSE file
Copyright (C) 2021 Michele Gandolfi, Michele Ceotto

Convert normal modes trajectory to Cartesian, given the rotation
matrix and masses.
Units must be consistent
Michele Gandolfi 2021
"""

import numpy as np
import sys
try: 
    from hessapprox import xyzfile
except Exception as e:
    #print( e)
    #print( 'Trying to temporarily add ".." to PYTHONPATH')
    sys.path.append( '..')
    try:
        from hessapprox import xyzfile
    except:
        quit( 'Cannot import module "xyzfile" for writing. Aborting..')
    else:
        pass
        #print( "Done")

assert len( sys.argv) == 6, 'You must pass 5 arguments: {0}  input  cnorm  mass  qrt  output' \
                        .format( sys.argv[0])
nmtraj, cnorm, mass, qrt, output = sys.argv[1:6]

b2a   = 0.529177210903
nmtraj= np.genfromtxt( nmtraj)
cnorm = np.genfromtxt( cnorm )
mass  = np.genfromtxt( mass  ).flatten()
qrt   = np.genfromtxt( qrt   ).flatten()

ncart = cnorm.shape[0]
cartmass = mass.repeat( 3)

cartraj = np.zeros( (nmtraj.shape[0], ncart))
for i, qvib in enumerate( nmtraj):
    q = np.hstack( ( qrt, qvib))
    cartraj[i,:]  = ( cnorm @ q ) / np.sqrt( cartmass) * b2a

cartraj = cartraj.reshape( -1, ncart//3, 3)

xyzfile.write_xyz( output, cartraj, ['O','H','H'])
