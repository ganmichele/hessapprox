#!/usr/bin/env python
"""
Locate optimal positions to compute the Hessian matrix given a trajectory
using the NGas method explained in
...

Michele Gandolfi 2021
"""

import numpy as np
import argparse 
import sys
from time import time

try:
    from hessapprox.Ngas import ngas
    from hessapprox import xyzfile
except Exception as e1:
    #print( 'Could not properly import hessapprox package. Received error with exception\n{0}'.format( e1))
    print( 'Assuming module is in ../hessapprox/. Trying to add to sys.path..')
    sys.path.append( '..')
    try:
        from hessapprox.Ngas import ngas
        from hessapprox import xyzfile
    except Exception as e:
        quit( 'Cannot import hessapprox package, Aborting..')
    #print( 'All packages imported successfully')

parser = argparse.ArgumentParser(
        description="""
                    Locate positions to compute the Hessian matrix, given a trajectory
                    using the NGas method explained in
                    """,
        formatter_class= argparse.ArgumentDefaultsHelpFormatter,
        fromfile_prefix_chars='+',
        epilog='Further arguments can be given in a file using +FILENAME')

# TRAJECTORY INPUT
parser.add_argument( 'trajectory' ,help='trajectory positions as a (multiple) space-separated values file (one position per line) or see the xyz flag')

# NEURAL GAS PARAMETERS
parser.add_argument( '-N', '--neurons'   ,default=100     ,help='number of neurons', type=int)
parser.add_argument( '-E', '--epochs'    ,default=30      ,help='Number of optimization epochs', type=int)
parser.add_argument( '-S', '--scale'     ,default='range' ,choices=['range','auto','var', 'cent']      ,help='scaling method to Neural Gas space')
parser.add_argument( '-I', '--init'      ,default='data'  ,choices=['data', 'ran', 'rdata', 'sphdata'] ,help='method to initialize neurons positions')
parser.add_argument( '--Lambda', nargs=2 ,default=None    ,help='custom lambda initial and final values')
parser.add_argument( '--Alpha' , nargs=2 ,default=None    ,help='custom alpha initial and final values')
parser.add_argument( '--neighbors'       ,default=None    ,help='custom number of neighbors*lambda for approximation', type=int)
parser.add_argument( '--dist'            ,default=None    ,help='custom distance function (NOT YET IMPLEMENTED)')
parser.add_argument( '--endon'           ,default='cent', choices=['traj', 'cent'] ,help='once the optimization is done, force the neurons to be on top of the trajectory ("traj", NOT RECOMMENDED) or at their center of mass ("cent")')
parser.add_argument( '--recenter_always' ,action='store_true' ,help='recenter neurons to their trajectory "center of mass" every epoch (NOT RECOMMENDED)')
parser.add_argument( '--moved_dist'      ,action='store_true' ,help='compute how much neurons have moved from their initial positions (in neurons space)')

# WHICH POINTS AND DEGREES OF FREEDOM YOU WANT TO CONSIDER
parser.add_argument( '-i', '--indexes'   ,default=['all'] ,nargs='+' ,help='indexes of trajectory to consider')
parser.add_argument( '-V','--var_indexes',default=['all'] ,nargs='+' ,help='indexes of variables (second axis) to consider')

# OUTPUTS & FORMATS
parser.add_argument( '-O', '--output', default='neurons', help='output file with neurons locations')
parser.add_argument( '-T', '--neurons_traj',action='store_true', help='outputs the neural gas approximate trajectory, instead of the neurons only')
parser.add_argument( '--xyz'         ,action='store_true' , help='input trajectory file is in xyz format')

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

# MORE PARSING
index    = _get_index_list( args.indexes)
varindex = _get_index_list( args.var_indexes)

if args.xyz:
    # file in Cartesian coordinates as xyz file
    atoms, xyz = xyzfile.read_xyz( args.trajectory, traj=True)
    n = len( atoms)*3
else:
    # file as space separated values (one geometry each row)
    xyz = np.genfromtxt( args.trajectory)
    n = xyz.shape[1]

not_varindex    = [e for e in np.arange( n) if e not in varindex]

xyz = xyz.reshape( -1, n)
if args.indexes[0] != 'all':
    xyz = xyz[index,:]

stationary_coor = xyz[:,not_varindex]

if args.var_indexes[0] != 'all':
    xyz = xyz[:,varindex]

# flatten the trajectory tensor to 2-D matrix so that each row is one geometry
kwargs = {}
if args.Lambda is not None:
    kwargs.update( {'l_i': float( args.Lambda[0]),
                    'l_f': float( args.Lambda[1])})
if args.Alpha is not None:
    kwargs.update( {'a_i': float( args.Alpha[0]),
                    'a_f': float( args.Alpha[1])})
if args.neighbors is not None:
    kwargs.update( {'nk' : int( args.neighbors)})

# NGAS
ngas = ngas( nn=args.neurons, **kwargs)

## if you have already trained the network and saved it you can load it with
## ngas.load_map( file), so you do not have to train it again from scratch
## the loaded map will not have the scaling parameter, you'd have to add them
## yourself

start = time()
ngas.train( xyz
           ,scale_method    = args.scale
           ,init_method     = args.init
           ,epochs          = args.epochs
           ,recenter_always = args.recenter_always
           ,endon           = args.endon
           ,moved_dist      = args.moved_dist
           )
end = time()
print( 'training time = {0}'.format( end-start))

if args.moved_dist:
    print( 'In optimization the neurons have moved: {0}'.format( ngas.moved_dist))

ngas_traj = ngas.scale_back() # if no arguments then neurons are scaled back
ngas_traj = np.hstack( (ngas_traj, stationary_coor[:ngas_traj.shape[0],:]))
if args.xyz:
    ngas_traj = ngas_traj.reshape( -1, len( atoms), 3) # Cartesian shape
    if args.output.split( '.')[-1:] != 'xyz':
        args.output += '.xyz'
    xyzfile.write_xyz( args.output, ngas_traj, atoms)
else:
    np.savetxt( args.output, ngas_traj)

bmus = ngas.get_BMUs( ngas.train_data).astype( int)
with open( 'relations_NGas.dat', 'w') as f:
    f.write( '# neurons-MDtraj relations. First index is 1 (not 0)\n')
    f.write( '# trajID    neuronID\n')
    for i, bmu in enumerate( bmus, start=1):
        f.write( '{0}\t\t{1}\n'.format( i, bmu+1))

if args.neurons_traj:
    # assign trajectory points to nearest neuron (Best Matching Unit or BMU):
    ngas_traj = ngas_traj[bmus,:] # expand over BMUs
    outputT   = args.output.rsplit('.', 1)[0]

    if args.xyz:
        ngas_traj = ngas_traj.reshape( -1, len( atoms), 3) # Cartesian shape
        xyzfile.write_xyz( outputT + '-traj.xyz', ngas_traj, atoms)
    else:
        np.savetxt( outputT + '-traj.dat', ngas_traj)

print( 'Neurons trajectory written in {0}'.format( args.output))
