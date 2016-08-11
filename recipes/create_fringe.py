#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 09:51
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
imcombine.py
"""

from common import open_with_meta

# from multiprocessing import Pool
# from multiprocessing import cpu_count
from common import combine_odi

import bottleneck as bn


def scale_to_bkg(data):
    mode = bn.nanmedian(data) * 3. - bn.nanmean(data) * 2
    data /= mode
    return data


def fringe_combine(data):
    return bn.nanmedian(data, axis=2) - 1


if __name__ == "__main__":
    import sys
    in_files = sys.argv[1:-1]
    out_file = sys.argv[-1]
    hdulist, exts, layout = open_with_meta(in_files[0])
    ret = combine_odi(in_files, fringe_combine, scale_to_bkg, exts=exts)
    for i, (ext, _) in enumerate(exts):
        hdulist[ext].data = ret[i]
    hdulist.writeto(out_file, clobber=True)
