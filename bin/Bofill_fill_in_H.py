#!/usr/bin/env python
"""
Hessian Approximation Methods. Full notice in LICENSE file
Copyright (C) 2021 Michele Gandolfi, Michele Ceotto

Use Hessian update scheme to fill in an incomplete list of
Hessian matrix along a MD simulation using Bofill method,
as explained in
H Wu, M Rahman, J Wang, U Louderaj, W L Hase, and Y Zhuang. J. Chem. Phys. 133, 074101 (2010) 

Michele Gandolfi 2021
"""

import numpy as np
import argparse 
import sys

try:
    from hessapprox.Ngas import ngas
    from hessapprox import xyzfile
except Exception as e1:
    print( 'Assuming module is in ../hessapprox/. Trying to add to sys.path..')
    sys.path.append( '..')
    try:
        from hessapprox.hessian_updt import fill_in_H
        from hessapprox import xyzfile
    except Exception as e:
        quit( 'Cannot import hessapprox package, Aborting..')

parser = argparse.ArgumentParser(
        description="""
                    Use Hessian update scheme to fill in an incomplete list of
                    Hessian matrix along a MD simulation using Bofill method, 
                    as explained in
                    H Wu, M Rahman, J Wang, U Louderaj, W L Hase, and Y Zhuang. J. Chem. Phys. 133, 074101 (2010) 
                    """,
        formatter_class= argparse.ArgumentDefaultsHelpFormatter,
        fromfile_prefix_chars='+',
        epilog='Further arguments can be given in a file using +FILENAME')

# TRAJECTORY INPUT
parser.add_argument( '-X', '--coor' ,help='trajectory positions as a (multiple) space-separated values file (one position per line)', required=True)
parser.add_argument( '-G', '--grad' ,help='gradients as a (multiple) space-separated values file (one gradient per line)', required=True)
parser.add_argument( '-H', '--hess' ,help='incomplete list of hessians as a (multiple) space-separated values file (one hessian per line)', required=True)

# HESSIAN UPDATE ARGUMENTS
parser.add_argument( '-N', '--nevals', default=10, help='number of hessians to predict, for each given', type=int)
parser.add_argument( '-L', '--lambd'             , help='Lambda value. Leave empty for Bofill optimal value', type=float)

# OUTPUT & FORMATS
parser.add_argument( '-O', '--output', default='hess_approx.dat', help='output file with DBq locations')
parser.add_argument( '--xyz'         ,action='store_true' , help='input files are in xyz format')

args = parser.parse_args()

if args.xyz:
    _, xyz = xyzfile.read_xyz( args.coor, traj=True)
    xyz = xyz.reshape( xyz.shape[0], -1)
else:
    xyz = np.genfromtxt( args.coor)

df    = xyz.shape[1]
grads = np.genfromtxt( args.grad)
hess  = np.genfromtxt( args.hess)
hessex= hess.copy()
hess  = hess.reshape( xyz.shape[0], hess.shape[1]//df, hess.shape[1]//df)
hess  = fill_in_H( xyz, grads, hess, step=args.nevals, lam=args.lambd, method=1).reshape( hess.shape[0], -1)
print( 'Difference between exact and approximate hess = {0}'.format( abs( hess - hessex).mean()))
tot = 0.0

np.savetxt( args.output, hess)
