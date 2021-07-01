#!/usr/bin/env python
"""
Locate positions to compute the Hessian matrix given a trajectory
with DBH method as explained in
...

Michele Gandolfi 2021
"""

import numpy as np
import argparse 
import sys

try:
    from hessapprox.DBH import DBH
    from hessapprox import xyzfile
except Exception as e1:
    #print( 'Could not properly import hessapprox package. Received error with exception\n{0}'.format( e1))
    print( 'Trying to add package to sys.path..')
    sys.path.append( '..')
    try:
        from hessapprox.DBH import DBH
        from hessapprox import xyzfile
    except Exception as e:
        quit( 'Cannot import hessapprox package, Aborting..')
    #print( 'All packages imported successfully')

parser = argparse.ArgumentParser(
        description="""
                    Locate positions to compute the Hessian matrix, given a trajectory with DBH method
                    """,
        formatter_class= argparse.ArgumentDefaultsHelpFormatter,
        fromfile_prefix_chars='+',
        epilog='Further arguments can be given in a file using +FILENAME')

# TRAJECTORY INPUT
parser.add_argument( 'trajectory' ,help='trajectory positions as a (multiple) space-separated values file (one position per line)')

# DBH GAS PARAMETER
parser.add_argument( '-R', '--rho', default=1.0, help='sphere radius', type=float)
parser.add_argument( '--dist'     , default=None,help='custom distance function (NOT YET IMPLEMENTED)')

# WHICH POINTS AND DEGREES OF FREEDOM YOU WANT TO CONSIDER
parser.add_argument( '-i', '--indexes'   ,default=['all'] ,nargs='+' ,help='indexes of trajectory to consider')
parser.add_argument( '-V','--var_indexes',default=['all'] ,nargs='+' ,help='indexes of variables (second axis) to consider')

# FILE NAMES, FORMATS AND TYPE OF OUTPUT
parser.add_argument( '-O', '--output', default='DBq_trajectory', help='output file with DBq locations')
parser.add_argument( '-T', '--DBq_traj',action='store_true', help='outputs the approximate trajectory, instead of the DBq only only')
parser.add_argument( '--xyz', action='store_true' , help='input trajectory file is in xyz format')

args = parser.parse_args()


def _get_index_list( ind):
    """
    Prepare index list from argument parser input
    """
    if ind[0]=='all':
        return {n for n in range( 10**5)}
    if ':' in ind[0]:
        splits = ind[0].split(':')
        if len( splits) == 3:
            low, step, high = splits
            low   = int( low)  if low  != '' else 0
            step  = int( step) if step != '' else 1
            high  = int( high)
            index = list( range( low, high, step))
        elif len( splits) == 2:
            low, high = splits
            low   = int( low)  if low  != '' else 0
            high  = int( high) #if high != '' else nmodes
            index = list( range( low, high))
        else:
            print( 'How many damn ":" are you using to define the indexes?!?!?!')
            quit()
    else:
        index = [int(i) for i in args.indexes]

    return index


# MORE PARSING OF OPTIONS
index    = _get_index_list( args.indexes)
varindex = _get_index_list( args.var_indexes)

if args.xyz:
    # file in Cartesian coordinates as xyz file
    atoms, xyz = xyzfile.read_xyz( args.trajectory, traj=True)
    n = len( atoms)*3
else:
    # file as space separated values (one geometry each row)
    #print( 'Assuming the trajectory is a file of space separated floats')
    xyz = np.genfromtxt( args.trajectory)
    n = xyz.shape[1]

not_varindex    = [e for e in np.arange( n) if e not in varindex]

xyz = xyz.reshape( -1, n)
if args.indexes[0] != 'all':
    xyz = xyz[index,:]

stationary_coor = xyz[:,not_varindex]

if args.var_indexes[0] != 'all':
    xyz = xyz[:,varindex]

# CREATE DBq 
DB = DBH()
DB.createDB( xyz, thr=args.rho)
dbq_traj = np.array( DB.DBq)

# WRITE RESULTS
dbq_traj = np.hstack( (dbq_traj, stationary_coor[:dbq_traj.shape[0],:]))
if args.xyz:
    dbq_traj = dbq_traj.reshape( -1, len( atoms), 3) # Cartesian shape
    if args.output.split( '.')[-1:] != 'xyz':
        args.output += '.xyz'
    xyzfile.write_xyz( args.output, dbq_traj, atoms)
else:
    np.savetxt( args.output, dbq_traj)

if args.DBq_traj:
    # assign trajectory points to nearest DBq (Best Matching Unit or BMU):
    #bmus = np.array( [e for l in DB.relations.values() for e in l])
    bmus = np.zeros( len( xyz)).astype( int)
    for k in DB.relations.keys():
        bmus[ DB.relations[k]] = k
    dbq_traj = dbq_traj[bmus,:] # expand over BMUs
    outputT   = args.output.rsplit('.', 1)[0]

    if args.xyz:
        dbq_traj = dbq_traj.reshape( -1, len( atoms), 3) # Cartesian shape
        xyzfile.write_xyz( outputT + '-traj.xyz', dbq_traj, atoms)
    else:
        np.savetxt( outputT + '-traj.dat', dbq_traj)

print( 'DBq trajectory written in {0} file'.format( args.output))

