#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-03 15:00
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
apply_segment_mask.py
"""


import numpy as np
from multiprocessing import Pool
from multiprocessing import cpu_count
from common import open_with_meta
import bottleneck as bn
import warnings
from apply_region_mask import apply_region_mask
from apply_region_mask import get_region_mask


def apply_segment_mask(data, seg, otaxy):
    print "working on OTA {0}".format(otaxy)
    data[seg > 0] = np.nan
    # also mask out the low value edges
    bkg = bn.nanmedian(data) * 3 - bn.nanmean(data) * 2
    std = bn.nanstd(data)
    with warnings.catch_warnings():
        warnings.filterwarnings(
                'ignore', 'invalid value encountered in less')
        data[np.isnan(data) | (data < bkg - 5 * std)] = np.nan
    return data


def mp_worker(args):
    data, seg, reg, hdr, otaxy = args
    data = apply_segment_mask(data, seg, otaxy)
    if reg is not None:
        data = apply_region_mask(data, reg, hdr, otaxy)
    return data


if __name__ == "__main__":
    import sys
    from astropy.io import fits

    reg_dir, in_file, seg_file, out_file = sys.argv[1:]
    hdulist, exts, layout = open_with_meta(in_file)
    seglist = fits.open(seg_file)
    print "focal plane layout: {0}".format(layout)

    # reg region file
    reg_file = get_region_mask(reg_dir, in_file)

    pool = Pool(cpu_count())
    ret = pool.map_async(
            mp_worker,
            [(hdulist[ext].data,
              seglist[ext].data,
              reg_file, hdulist[ext].header,
              otaxy) for ext, otaxy in exts]
            ).get(9999999)
    for i, (ext, _) in enumerate(exts):
        hdulist[ext].data = ret[i]
    hdulist.writeto(out_file, clobber=True)
