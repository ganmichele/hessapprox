#!/usr/bin/env python
"""
DBH method for Hessian matrix approximation along MD trajectories

Michele Gandolfi 2021
"""

import numpy as np
import sys

def set_set_distance( x, y, p=2):
    n, m = x.shape[0], y.shape[0]
    D = np.zeros( (n,m))
    for i, a in enumerate( x):
        D[i,:] = set_point_distance( y, a, p=p)
        #for j, b in enumerate( y):
        #    D[i,j] = distance( a, b)
    return D

try:
    from scipy.spatial import distance_matrix as set_set_distance
except Exception as e:
    print( 'Cannot import scipy fast distance calculation')
    print( 'Use builtin slow distance calculator')


def set_point_distance( x, b, p=2):
    n = len( x)
    D = np.zeros( n)
    for i, a in enumerate( x):
        D[i] = distance( a, b, p=p)
    return D

def distance( a, b, p=2):
    a, b = np.array( a), np.array( b)
    if p==2:
        res = np.sqrt( (a-b) @ (a-b))
    elif p==1:
        res = abs( a - b).sum()   
    elif p==-1:
        res = abs( a - b).max()
    return res


class DBH():
    """
    Hessian DataBase (DBH) class
    
    METHODS:
    self.createDB( data, thr)
    self.updateDB( data)
    self.readDB( x)
    self._inDB( x)
    """
    def __init__( self, Pes=None):
        """
        Initialize an empty database
        """
        self.DBq = []
        self.DBH = []
        self.relations = {}
        if Pes is not None:
            self.Pes = Pes

    def _inDB( self, x):
        """
        Check if x is within self.thr from any point inside the DBq
        """
        if 'scipy' in sys.modules and len( self.DBq) > 0:
            D = set_set_distance( np.array(self.DBq), x.reshape(1,-1), p=self.p).flatten()
        else:
            D = set_point_distance( self.DBq, x, p=self.p)

        if np.any( D < self.thr):
            return True
        else:
            return False

    def updateDB( self, x, H=None, **pes_args):
        """
        Update DBq and, possibly DBH
        """
        self.DBq.append( x)
        if hasattr( self, 'Pes') and H == 'compute':
            H = self.Pes.hessian( x, **pes_args)
        elif H is not None:
            self.DBH.append( H)
            
    def createDB( self, data, thr=5.0, dist_measure='th', dohess=None, **pes_args):
        """
        Fill in the DB all the geometries so that none is further than "thr"
        from all points in the database
        """
        self.data_input = data
        self.thr = thr
        self.dist_measure = dist_measure
        if self.dist_measure[:2] == 'cb':
            self.p = 1
        elif self.dist_measure[:2] == 'eu':
            self.p = 2
        elif self.dist_measure[:2] == 'th':
            self.p = -1 if 'scipy' not in sys.modules else np.infty
        else: # assuming theta (minkowski p=infty distance)
            self.p = -1 if 'scipy' not in sys.modules else np.infty

        for x in data:
            if not self._inDB( x):
                self.updateDB( x, H=dohess, **pes_args)

        self.update_relations()

    def update_relations( self):
        #dist_q_DBq = set_set_distance( np.array( self.DBq), 
        #                               np.array( self.data_input))
        self.relations = {k: [] for k in range( len( self.DBq))}

        if 'scipy' in sys.modules:
            D = set_set_distance( np.array(self.DBq), self.data_input, p=self.p)
            K = D.argmin(0)
            for i, j in enumerate( K):
                self.relations[j].append( i)

        else:
            for i, q in enumerate( self.data_input):
                dbID = np.argmin(
                    set_point_distance( np.array( self.DBq), 
                                        q, 
                                        #measure=self.dist_measure)
                                        #measure='eucl'
                                        p=self.p)
                    )
                self.relations[dbID].append( i)


if __name__ == '__main__':
    print( 'Nothing to do')
    pass
