#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-30 23:04
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
apply_ota_mask.py
"""


import os
import glob
import re
import numpy as np
from multiprocessing import Pool
from multiprocessing import cpu_count
from common import open_with_meta


def apply_bpmask(data, otaxy, layout):
    if layout == 'odi56':
        bpmdir = '/home/ma/Codes/podi/trunk/.bpm/odi_5x6'
    elif layout == 'podi':
        bpmdir = '/home/ma/Codes/podi/trunk/.bpm/podi/'
    else:
        raise RuntimeError('focal plane layout {0} not recognized'
                           .format(layout))
    bpm_file = os.path.join(bpmdir, 'bpm_xy{0}.reg'.format(otaxy))
    bpm = []
    with open(bpm_file, 'r') as fo:
        for ln in fo.readlines():
            rect = re.match(r'box\(([0-9+-., ]+)\)', ln.strip())
            if rect is not None:
                rect = map(float, rect.group(1).split(','))
                # print "box from bpm: {0}".format(rect)
                bpm.append((
                    rect[0] - rect[2] * 0.5,
                    rect[0] + rect[2] * 0.5,
                    rect[1] - rect[3] * 0.5,
                    rect[1] + rect[3] * 0.5))
            else:
                continue
    for box in bpm:
        l, r, b, t = [int(v + 0.5) for v in box]
        data[b:t, l:r] = np.nan
    return data


def get_mask_ota(in_file, flag_file):
    otas = []
    with open(flag_file, 'r') as fo:
        for ln in fo.readlines():
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            ln = ln.strip().split()
            if len(ln) < 2:
                continue
            else:
                if ln == '*':
                    otas.extend(map(int, ln[1:]))
                else:
                    match = [f for f in glob.glob(ln[0])
                             if os.path.samefile(f, in_file)]
                    if len(match) == 1:
                        if ln[1] == '*':
                            otas.append(-1)
                        else:
                            otas.extend(map(int, ln[1:]))
                    elif len(match) == 0:
                        continue
                    else:
                        raise RuntimeError('unable to get OTA mask for {0}'
                                           .format(in_file))
    return otas


if __name__ == "__main__":

    import sys

    in_file, flag_file, out_file = sys.argv[1:]
    otas = get_mask_ota(in_file, flag_file)
    if -1 in otas:
        print "mask ota {0} for {1}, discard".format(otas, in_file)
        sys.exit(0)
    print "mask ota {0} for {1}".format(otas, in_file)

    hdulist, exts, layout = open_with_meta(in_file)
    print "focal plane layout: {0}".format(layout)

    def mp_worker(args):
        data, otaxy = args
        print "working on OTA {0}".format(otaxy)
        if otaxy in otas:
            data[:, :] = np.nan
        else:
            data = apply_bpmask(data, otaxy, layout)
        return data

    pool = Pool(cpu_count())
    ret = pool.map_async(
            mp_worker,
            [(hdulist[i].data, otaxy) for i, otaxy in exts]
            ).get(9999999)
    for ii, (i, _) in enumerate(exts):
        hdulist[i].data = ret[ii]
    hdulist.writeto(out_file, clobber=True)
