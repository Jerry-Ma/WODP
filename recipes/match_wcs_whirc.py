#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-02 02:52
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
match_wcs_whirc.py

Match sextractor catalogs to a common wcs

group the sources, such that each group contains one real source

matrix
    c    c    c    c    c
s
s   x    x    x
s   x    x         x    x
s
s

as this quantity is the least

std(x[s, c+1] - x[s, c])

"""

import numpy as np
import contextlib


@contextlib.contextmanager
def printoptions(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    yield
    np.set_printoptions(**original)


if __name__ == "__main__":
    import sys
    from astropy.table import Table
    in_files = sys.argv[1:-1]
    out_file = sys.argv[-1]
    ncat = len(in_files)
    nslot = 10
    mat = np.ones((nslot, ncat, 3), dtype='d') * np.nan
    for i, f in enumerate(in_files):
        tbl = Table.read(f, format='ascii.sextractor')
        tbl.sort("FLUX_AUTO")
        tbl.reverse()
        mat[:, i, 0] = tbl['MAG_AUTO'][:nslot]
        mat[:, i, 1] = tbl['X_IMAGE'][:nslot]
        mat[:, i, 2] = tbl['Y_IMAGE'][:nslot]
    np.set_printoptions(
            formatter={'float': '{:10.3f}'.format},
            linewidth=np.inf)
    mm = np.ones((nslot, ncat))
    ii = np.repeat(np.arange(nslot), ncat).reshape((nslot, ncat))
    jj = np.repeat(np.arange(ncat), nslot).reshape((ncat, nslot)).T

    # cost = squared sum over all two selected j1, j2
    # for for each element, three states:
    #   masked, unmasked, swapped with other
    # loop over the sates, and get the min final lost
    jlist = range(0, ncat)
    j = 0

    def get_dcost(j1, j2, sv, v):
        ii1 = np.copy(ii[:, j1])
        ii2 = np.copy(ii[:, j2])
        m1 = np.copy(mm[:, j1][ii1])
        m2 = np.copy(mm[:, j2][ii2])
        # do the swap or mask
        if v == -1:
            m1[sv] = np.nan
        elif np.isnan(m1[sv]) and np.isnan(m1[v]):
            return np.inf
        else:
            tmp = ii1[sv]
            ii1[sv] = ii1[v]
            ii1[v] = tmp
        x1 = mat[:, j1, 1][ii1]
        x2 = mat[:, j2, 1][ii2]
        y1 = mat[:, j1, 2][ii1]
        y2 = mat[:, j2, 2][ii2]
        cost = np.nanvar(
                np.hypot((x1 - x2) * m1, (y1 - y2) * m2))
        return cost

    def get_cost(jlist):
        cost = 0
        for cj in jlist:
            subjlist = [j for j in jlist if j != cj]
            dcost = np.inf
            oper = (0, 0)
            for sv in range(nslot):
                for v in range(-1, nslot):
                    _dcost = sum([get_dcost(cj, j, sv, v) for j in subjlist])
                    if _dcost < dcost:
                        dcost = _dcost
                        oper = (sv, v)
            if oper == (0, 0):
                print "stable"
                pass
            else:
                sv, v = oper
                if v == -1:
                    mm[sv, cj] = np.nan
                else:
                    tmp = ii[sv, cj]
                    ii[sv, cj] = ii[v, cj]
                    ii[v, cj] = tmp
            if len(subjlist) == 1:
                cost += dcost
            else:
                cost += dcost + get_cost(subjlist)
        print mm
        print ii
        if cost = 0:
        return cost
    get_cost(range(ncat))
    print (mat[:, :, 0] * mm)[ii, jj]
