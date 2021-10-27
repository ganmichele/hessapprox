#!/usr/bin/env python
"""
Hessian Approximation Methods. Full notice in LICENSE file
Copyright (C) 2021 Michele Gandolfi, Michele Ceotto

DBH method for Hessian matrix approximation along MD trajectories.
For reference see the following articles:
Conte R, Gabas F, Botti G, Zhuang Y, Ceotto M J. Chem. Phys. 150, 244118 (2019)
Gandolfi M, Ceotto M, J. Chem. Theory Comput. (2021) (minor modification of the original method)

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
    self.createDB
    self.updateDB
    self.update_relations
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

    def createDB(  self, data
                  ,thr=5.0
                  ,dist_measure='th'
                  ,dohess=None
                  ,bestrel=True
                  ,**pes_args):
        """
        Create the DB given the data array and threshold thr.
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

        for xid, x in enumerate( data):
            self.updateDB( x, xid, H=dohess, **pes_args)

        if bestrel:
            self.update_relations()


    def _inDB( self, x):
        """
        Check if x is within self.thr from any point inside the DBq
        """
        if 'scipy' in sys.modules and len( self.DBq) > 0:
            D = set_set_distance( np.array(self.DBq), x.reshape(1,-1), p=self.p).flatten()
        else:
            D = set_point_distance( self.DBq, x, p=self.p)

        if np.any( D < self.thr):
            self._who = np.argmin( D)
            return True
        else:
            return False


    def updateDB( self, x, xid, H=None, **pes_args):
        """
        Update DBq and, possibly DBH if H is passed
        xid is the index for the trajectory entry x.
        If H=='compute' and the Pes attribute was passed, the
        Hessian matrix is computed on-the-fly (assuming that Pes
        computes the hessian with the 'hessian' method).
        """
        if not self._inDB( x):
            self.relations[len( self.DBq)] = [xid,]
            self.DBq.append( x)
            if hasattr( self, 'Pes') and H == 'compute':
                H = self.Pes.hessian( x, **pes_args)
            elif H is not None:
                self.DBH.append( H)
        else:
            self.relations[self._who].append( xid)
            

    def update_relations( self, alldata=None):
        """
        Wipe down current relations and create new ones.
        To be used after the last update.
        Optional argument alldata contains the configurations of you trajectory.
        If alldata is not passed, self.data_input will be used.
        """
        if alldata:
            self.data_input = alldata

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
