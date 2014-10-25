# -*- coding: utf-8 -*-
"""
Created on Tue May 15 10:47:45 2012

@author: nelson
"""
import pylab as P


def t3dinp(fname='hitpopn.01.t3d'):
    fid = open(fname)
    # get header
    # skip first line
    fid.readline()
    # return first three of next four values
    np, ne, nt = P.array(fid.readline().split()[:-1], dtype=P.int32)
#    print np, ne, nt

    # load rest of data
    arr = P.loadtxt(fid, skiprows=1)

    r = arr[:np, 1:3]
    pt = 1000 * arr[:np, 4] + arr[:np, 5]
    # Python starts indexes at 0
    tris = (arr[np:, 1:4] - 1).astype(P.int32)

    return r, tris, pt