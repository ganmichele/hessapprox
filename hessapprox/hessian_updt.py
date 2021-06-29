#!/usr/bin/env python
"""
Update Hessian matrix according to Bofill's method.
Michele Gandolfi 2021
as in
Gandolfi M, Ceotto M. title . paper . 2021
"""

import numpy as np


def HessUpd( dq, dg, H, lam=None):#, version=2):
    """
    Update the Hessian matrix to H(q2), given positions and
    gradient at positions q1 and q2.
    INPUT:
    dq = q2 - q1    coordinate difference
    dg = g2 - g1    gradient difference
    H               Hessian at q1
    lam             lambda factor, if None, default (Bofill) value 
                    is used
    """
    Hout = H.copy()
    R = 2.0 * ( dg - H @ dq)
    RoR  = np.outer( R , R )
    dqoR = np.outer( dq, R )
    dqodq= np.outer( dq, dq)
    RR   = np.dot  ( R , R )
    Rdq  = np.dot  ( R , dq)
    dqdq = np.dot  ( dq, dq)

    if lam is None:
        lam = 1.0 - (Rdq**2 / ( RR * dqdq) )

    Hout += (1-lam) * RoR / Rdq
    Hout += lam * ( dqoR + dqoR.T ) / dqdq
    Hout -= lam * Rdq / dqdq**2 * dqodq
    return Hout


def fill_in_H( qarr, garr, Harr, step=10, lam=None, method=1):
    """
    fill in the missing hessians in Harr using Bofill hessian updates
    method=1 (default) means that the last approximated hessian is used
    to extrapolate the next one (so q2 is the successive step of q1).
    If method == 2 then q2 and q1 might not be successive steps.
    """
    assert qarr.shape == garr.shape, 'positions and gradients must have same shape'
    assert qarr.shape[0] == Harr.shape[0], 'hessian array must contain same number of entries as qarr'
    assert qarr.shape[1] == Harr.shape[1], 'invalid hessian matrix shape'
    
    Harr_apprx = Harr.copy()
    if method not in {1,2}: method=1

    if method==1:
    # Use last approximated hessian to extrapolate next (default)
        dqarr = np.diff( qarr, axis=0)
        dgarr = np.diff( garr, axis=0)
        for i, (dq, dg) in enumerate( zip( dqarr, dgarr), start=1):
            if i % step != 0:
                Harr_apprx[i] = HessUpd( dq, dg, Harr_apprx[i-1], lam=lam)
    elif method==2:
    # Use last exact hessian to extrapolate next
        ref = 0
        for i in range( 1, len( qarr)):
            dq = qarr[i] - qarr[ref]
            dg = garr[i] - garr[ref]
            if i % step != 0 or i == 0:
                Harr_apprx[i] = HessUpd( dq, dg, Harr_apprx[ref], lam=None)
            else:
                ref = i

    return Harr_apprx

