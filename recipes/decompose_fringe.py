#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 11:55
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
decompose_fringe.py
"""

from common import open_with_meta
import numpy as np
from pyjerry.instrument import wiyn
import itertools
from scipy import interpolate


def decompose(data):
    '''decompose the pattern to HF and LF'''
    wl = wiyn.WIYNLayout(binning=1.)
    # bin each cell by nbin
    nbin = 8
    bw = wl.CW / nbin
    bh = wl.CH / nbin
    # bs = np.hypot(bw, bh) * 0.5
    samp_v = np.empty((wl.NCY * nbin, wl.NCX * nbin))
    samp_i = []
    samp_j = []
    for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        (cl, cr), (cb, ct) = wl.get_cell_rect(0, 0, cj, ci)
        for bj, bi in itertools.product(range(nbin), repeat=2):
            bl = cl + bj * bw
            br = bl + bw
            bb = cb + bi * bh
            bt = bb + bh
            bx = 0.5 * (bl + br)
            by = 0.5 * (bb + bt)
            v = np.nanmedian(data[bb:bt, bl:br])
            samp_v[ci * nbin + bi, cj * nbin + bj] = v
            if cj == ci and bj == bi:
                samp_j.append(bx)
                samp_i.append(by)
    samp_i, samp_j = np.meshgrid(np.array(samp_i),
                                 np.array(samp_j),
                                 indexing='ij')
    # smooth the anchor array before spline
    # box_2D_kernel = Box2DKernel(3)
    # samp_v = convolve(samp_v, box_2D_kernel, boundary='extend')
    # samp_v = samp_v - np.nanmedian(samp_v)
    m = ~np.isnan(samp_v)
    kx = ky = 3
    print "size of sampling array: {0} {1}".format(samp_v[m].size, samp_v.size)
    if len(samp_v[m]) < (kx + 1) * (ky + 1):
        raise RuntimeError("unable to create LF map")
    else:
        # do spline interpolation
        spline = interpolate.bisplrep(
                samp_i[m], samp_j[m], samp_v[m], kx=kx, ky=ky)
        ii, jj = np.mgrid[0:data.shape[0], 0:data.shape[1]]
        _lf = interpolate.bisplev(ii[:, 0], jj[0, :], spline).reshape(
                data.shape)
        # mask bad cell off on lf
        lf = np.ones_like(_lf) * np.nan
        for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
            (cl, cr), (cb, ct) = wl.get_cell_rect(0, 0, cj, ci)
            if np.any(~np.isnan(data[cb:ct, cl:cr])):
                lf[cb:ct, cl:cr] = _lf[cb:ct, cl:cr]
    hf = data - lf
    return hf, lf


if __name__ == "__main__":
    import sys
    in_file, hf_file, lf_file = sys.argv[1:]
    hdulist, exts, layout = open_with_meta(in_file)
    hfdata = []
    lfdata = []
    for ext, otaxy in exts:
        print "working on OTA {0}".format(otaxy)
        hf, lf = decompose(hdulist[ext].data)
        # make use of the first images as template
        hfdata.append(hf)
        lfdata.append(lf)
    for i, (ext, _) in enumerate(exts):
        hdulist[ext].data = hfdata[i]
    hdulist.writeto(hf_file, clobber=True)
    for i, (ext, _) in enumerate(exts):
        hdulist[ext].data = lfdata[i]
    hdulist.writeto(lf_file, clobber=True)
    hdulist.close()
