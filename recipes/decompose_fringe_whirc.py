#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 11:55
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
decompose_fringe.py
"""

import numpy as np
import itertools
from scipy import interpolate
from astropy.convolution import Box2DKernel
from astropy.convolution import convolve


def decompose(data):
    '''decompose the pattern to HF and LF'''
    nbin = 8
    bw = data.shape[1] / nbin
    bh = data.shape[0] / nbin
    # bs = np.hypot(bw, bh) * 0.5
    samp_v = np.empty((nbin, nbin))
    samp_i = []
    samp_j = []
    for bj, bi in itertools.product(range(nbin), repeat=2):
        bl = bj * bw
        br = bl + bw
        bb = bi * bh
        bt = bb + bh
        bx = 0.5 * (bl + br)
        by = 0.5 * (bb + bt)
        v = np.nanmedian(data[bb:bt, bl:br])
        samp_v[bi, bj] = v
        if bj == bi:
            samp_j.append(bx)
            samp_i.append(by)
    samp_i, samp_j = np.meshgrid(np.array(samp_i),
                                 np.array(samp_j),
                                 indexing='ij')
    # smooth the anchor array before spline
    box_2D_kernel = Box2DKernel(3)
    samp_v = convolve(samp_v, box_2D_kernel, boundary='extend')
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
        lf = interpolate.bisplev(ii[:, 0], jj[0, :], spline).reshape(
                data.shape)
    hf = data - lf
    return hf, lf


if __name__ == "__main__":
    import sys
    from astropy.io import fits
    in_file, hf_file, lf_file = sys.argv[1:]
    hdulist = fits.open(in_file)
    exts = [(0, 'WHIRC')]
    hfdata = []
    lfdata = []
    for ext, extname in exts:
        print "working on frame {0}".format(extname)
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
