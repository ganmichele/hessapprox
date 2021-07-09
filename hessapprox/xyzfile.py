#!/usr/bin/env python
"""
Hessian Approximation Methods. Full notice in LICENSE file
Copyright (C) 2021 Michele Gandolfi, Michele Ceotto

Read write and interpret coordinates files in xyz format

Michele Gandolfi 2021
"""

import numpy as np
import re

def get_mass( atoms):
    """
    Get list of atoms' atomic weights
    """
    masses = []
    #for a in atoms:
    #    masses.append( mendeleev.element( a).atomic_weight)
    #return masses
    massdic = { 'H': 1837.1500,
                'D': 3674.3000,
                'C': 21874.660,
                'N': 25726.553,
                'O': 29156.9600,
                'S': 58450.9200}
    for a in atoms:
        masses.append( massdic[ a])
    return masses

def read_xyz( fname, traj=False, vel=False, until=None):
    """
    Read xyz file format and return useful content
    """
    regex_natoms = r'\d+'
    regex_split_coor = r'\s+'
    with open( fname, 'r') as f:
        s = f.readline() # number of atoms
        natoms = int( re.search( regex_natoms, s).group())
        atoms = []
        #xyz = np.zeros( (natoms, 3))
        xyz = []
        velocity = []
        f.readline() # whatever caption line 
        if until is not None:
            until *= ( natoms + 2)
        else:
            until = np.infty # READ EVERYTHING
        for i, l in enumerate( f):
            lclean = l.strip()
            if i+1 > natoms and not traj:
                if lclean != '':
                    print( 'Warning, detected more lines then declared atoms')
                    break
                else:
                    break
            else:
                if i>1 and i % (natoms+2) in {natoms,natoms+1}:
                    if i >= until:
                        break
                    else:
                        continue
            values = re.split( regex_split_coor, lclean.strip())
            if i < natoms:
                atoms.append( values[0])
            xyz += [ float( v) for v in values[1:4]]
            if vel:
                velocity += [float( v) for v in values[4:7]]
            #xyz[i,0] = float( values[1])
            #xyz[i,1] = float( values[2])
            #xyz[i,2] = float( values[3])
    xyz = np.array( xyz)
    xyz = xyz.reshape( -1, 3) if not traj else xyz.reshape( -1, natoms, 3)
    assert natoms == len( atoms), 'Something odd in this xyz file'
    if vel:
        velocity = np.array( velocity)
        velocity = velocity.reshape( -1, 3) if not traj \
                                            else velocity.reshape( -1, natoms, 3)
        return atoms, xyz, velocity
    return atoms, xyz

def write_xyz( fname, xyz, atoms, append=False):
    """
    Write xyz single point or trajectory to file
    """
    how_write = 'a' if append else 'w'
    with open( fname, how_write) as f:
        if len( xyz.shape) == 3:
            print( 'Writing trajectory to file {0}..'.format( fname))
            for p in xyz:
                f.write( '{0}\n'.format( p.shape[0])) # number of atoms
                f.write( '# atom x y z\n')            # blanck comment line
                for i, a in enumerate( p):
                    f.write( '{0} {1} {2} {3}\n'.format( atoms[i], a[0], a[1], a[2]))
        elif len( xyz.shape) == 2:
            f.write( '{0}\n'.format( xyz.shape[0])) # number of atoms
            f.write( '# atom x y z\n')            # blanck comment line
            for i, a in enumerate( xyz):
                f.write( '{0} {1} {2} {3}\n'.format( atoms[i], a[0], a[1], a[2]))
        else:
            print( 'Coordinates not understood.\nYou should pass a 2-D or 3-D numpy array as xyz coordinates')
            

if __name__ == '__main__':
    print( 'Nothing to do')
    pass
