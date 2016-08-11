#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 12:55
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
remove_fringe.py
"""

import os
import re
import numpy as np
import bottleneck
from astropy.stats import sigma_clipped_stats


def de_fringe(data, pattern, hf, otaxy):
    # skip all nan extension
    if np.all(np.isnan(data)):
        print "skip OTA {0}".format(otaxy)
        return data
    # measure amplitude
    # 'fringevectors__odi_i__OTA11.reg'
    fringe_vector_dir = '/home/ma/Codes/podi/trunk/.fringevectors/odi_5x6/'
    regfile = os.path.join(
            fringe_vector_dir,
            'fringevectors__odi_z__OTA{0:d}.reg'.format(otaxy))
    re_vector = r'.+vector\(\s*' \
                r'([+-]?\d+(?:\.\d+)?)\s*,\s*' \
                r'([+-]?\d+(?:\.\d+)?)\s*,\s*' \
                r'([+-]?\d+(?:\.\d+)?)\s*,\s*' \
                r'([+-]?\d+(?:\.\d+)?)' \
                r'\s*\).*'
    boxsize = 5
    samp = []
    print "read in vector file: {0}".format(regfile)
    with open(regfile, 'r') as fo:
        for ln in fo.readlines():
            matched = re.match(re_vector, ln)
            if matched is not None:
                # print "vector from region: {0}".format(ln)
                x1 = int(float(matched.group(1)) + 0.5)
                y1 = int(float(matched.group(2)) + 0.5)
                length = float(matched.group(3))
                angle = float(matched.group(4)) * np.pi / 180.
                x2 = int(x1 + length * np.cos(angle) + 0.5)
                y2 = int(y1 + length * np.sin(angle) + 0.5)
                d_d = bottleneck.nanmedian(
                        data[y1 - boxsize:y1 + boxsize,
                             x1 - boxsize:x1 + boxsize])
                d_l = bottleneck.nanmedian(
                        data[y2 - boxsize:y2 + boxsize,
                             x2 - boxsize:x2 + boxsize])
                p_d = bottleneck.nanmedian(
                        hf[y1 - boxsize:y1 + boxsize,
                           x1 - boxsize:x1 + boxsize])
                p_l = bottleneck.nanmedian(
                        hf[y2 - boxsize:y2 + boxsize,
                           x2 - boxsize:x2 + boxsize])
                d_a = d_l - d_d
                p_a = p_l - p_d
                if (d_a != 0) and (p_a != 0):
                    scale = d_a / p_a
                    samp.append([d_d, d_l, d_a, p_d, p_l, p_a, scale])
    samp = np.array(samp)
    scale = samp[:, -1]
    _, med, std = sigma_clipped_stats(scale, sigma=3.0)
    print "scaling factor: {0} +/- {1}".format(med, std)
    data -= pattern * med
    return data


if __name__ == "__main__":

    from astropy.io import fits
    import sys
    from common import open_with_meta

    in_file, pattern_file, hf_file, out_file = sys.argv[1:]

    hdulist, exts, _ = open_with_meta(in_file)
    pattern = fits.open(pattern_file)
    hf = fits.open(hf_file)

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
        hdulist[ext].data = de_fringe(
                hdulist[ext].data,
                pattern[ext].data,
                hf[ext].data,
                otaxy
                )
    hdulist.writeto(out_file, clobber=True)
    hdulist.close()
    pattern.close()
    hf.close()
