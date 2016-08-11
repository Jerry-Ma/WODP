#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-02 23:53
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
create_illum.py
"""


from common import open_with_meta
from common import combine_odi
from pyjerry.instrument import wiyn
# from astropy.convolution import convolve, Box2DKernel
from multiprocessing import Pool
from multiprocessing import cpu_count
import bottleneck as bn
import itertools
import numpy as np
from astropy.stats import sigma_clip


def scale_to_bkg(data):
    mode = bn.nanmedian(data) * 3. - bn.nanmean(data) * 2
    data /= mode
    return data


def illum_combine(data):
    return bn.nanmedian(data, axis=2)


def smooth_tile(data, otaxy, width):

    print "smoothing OTA {0}".format(otaxy)
    min_count = width * 0.1
    wl = wiyn.WIYNLayout(binning=1)
    # 2sigma clip
    for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        # print "smoothing OTA {0} cell {1}{2}".format(otaxy, cj, ci)
        (l, r), (b, t) = wl.get_cell_rect(0, 0, cj, ci)
        cell = data[b:t, l:r]
        if np.all(np.isnan(cell)):
            print "cell {0}{1} skip due to all NAN".format(cj, ci)
            continue
        cell = sigma_clip(cell, sigma=2,
                          cenfunc=bn.nanmedian, stdfunc=bn.nanstd)
        data[b:t, l:r] = cell
    _data = bn.move_mean(data, width, min_count=min_count, axis=1)
    _data = bn.move_mean(_data, width, min_count=min_count, axis=0)
    for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        # print "smoothing OTA {0} cell {1}{2}".format(otaxy, cj, ci)
        (l, r), (b, t) = wl.get_cell_rect(0, 0, cj, ci)
        if np.all(np.isnan(data[b:t, l:r])):
            print "cell {0}{1} skip due to all NAN".format(cj, ci)
            continue
        else:
            data[b:t, l:r] = _data[b:t, l:r]
    # cell = data[b:t, l:r]
    # data[b:t, l:r] = cell
    # return data
    return data


def mp_worker(args):
    data, otaxy, width = args
    return smooth_tile(data, otaxy, width)


if __name__ == "__main__":
    import sys
    in_files = sys.argv[1:-1]
    out_file = sys.argv[-1]
    hdulist, exts, layout = open_with_meta(in_files[0])
    ret = combine_odi(in_files, illum_combine, scale_to_bkg, exts=exts,
                      jobs=1)
    for i, (ext, _) in enumerate(exts):
        hdulist[ext].data = ret[i]
    hdulist.writeto(out_file.replace(".fits", '.check'), clobber=True)
    # do the smooth
    mp = Pool(cpu_count())

    width = 1000
    # kernel = Box2DKernel(width)
    ret = []
    for ext, otaxy in exts:
        ret.append(smooth_tile(hdulist[ext].data, otaxy, width))
    # ret = mp.imap_async(mp_worker, ((hdulist[ext].data, otaxy, width)
    #                                 for ext, otaxy in exts)).get(9999999)
    for i, (ext, _) in enumerate(exts):
        hdulist[ext].data = ret[i]
    hdulist.writeto(out_file, clobber=True)
