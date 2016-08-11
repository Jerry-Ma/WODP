#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 12:55
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
correct_illum.py
"""

import numpy as np


def correct(data, illum, otaxy):
    # skip all nan extension
    if np.all(np.isnan(data)):
        print "skip OTA {0}".format(otaxy)
        return data
    data /= illum
    return data


if __name__ == "__main__":

    from astropy.io import fits
    import sys
    from common import open_with_meta

    in_file, illum_file, out_file = sys.argv[1:]

    hdulist, exts, _ = open_with_meta(in_file)
    illum = fits.open(illum_file)

    # from multiprocessing import Pool
    # from multiprocessing import cpu_count

    # def mp_worker(args):
    #     data, pattern, hf, otaxy = args
    #     print "working on OTA {0}".format(otaxy)
    #     return de_fringe(data, pattern, hf, otaxy)

    # pool = Pool(cpu_count())
    # ret = pool.map_async(
    #         mp_worker,
    #         [(hdulist[ext].data, pattern[ext].data, hf[ext].data, otaxy)
    #          for ext, otaxy in exts]
    #         ).get(9999999)
    # for i, (ext, _) in enumerate(exts):
    #     hdulist[ext] = ret[i]

    for ext, otaxy in exts:
        hdulist[ext].data = correct(hdulist[ext].data, illum[ext].data, otaxy)
    hdulist.writeto(out_file, clobber=True)
    hdulist.close()
    illum.close()
