#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-03 16:02
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
apply_gsc_mask.py
"""


import os
import glob
import re
import numpy as np
from multiprocessing import Pool
from multiprocessing import cpu_count
from common import open_with_meta
import pyregion
from astropy.io import fits


def apply_region_mask(data, reg_file, hdr, otaxy):
    print "apply mask region {0}".format(reg_file)
    # create a fresh wcs from keys
    hdu = fits.ImageHDU(data=data)
    # add basic wcs
    for key in ['EQUINOX', 'CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2',
                'CRVAL1', 'CRVAL2', 'CUNIT1', 'CUNIT2', 'CD1_1', 'CD2_1',
                'CD1_2', 'CD2_2']:
        hdu.header[key] = hdr[key]
    mask = pyregion.open(reg_file).get_mask(hdu=hdu)
    data[mask] = np.nan
    return data


def get_region_mask(reg_dir, in_file):
    # looking for reg file
    regex = r'20\d{6}T\d{6}\.\d'
    id_ = re.search(regex, in_file).group(0)
    reg_files = [f for f in glob.glob(os.path.join(reg_dir, '*.reg'))
                 if id_ in f]
    # print id_, reg_files
    if len(reg_files) > 1:
        raise RuntimeError("multiple mask regions found for {0}"
                           .format(in_file))
    elif len(reg_files) == 0:
        return None
    else:
        return reg_files[0]


def mp_worker(args):
    apply_region_mask(*args)


if __name__ == "__main__":

    import sys

    reg_dir, in_file, out_file = sys.argv[1:]
    hdulist, exts, layout = open_with_meta(in_file)
    print "focal plane layout: {0}".format(layout)
    reg_file = get_region_mask(reg_dir, in_file)
    if reg_file is None:
        print "no mask region found for {0}".format(in_file)
        hdulist.writeto(out_file, clobber=True)
    else:
        print "object mask region: {0}".format(reg_file)
        pool = Pool(cpu_count())
        ret = pool.map_async(
                mp_worker,
                [(hdulist[ext], reg_file, otaxy) for ext, otaxy in exts]
                ).get(9999999)
        for i, (ext, _) in enumerate(exts):
            hdulist[ext].data = ret[i]
        hdulist.writeto(out_file, clobber=True)
