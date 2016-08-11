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


def apply_gsc(hdu, gsc, layout):



if __name__ == "__main__":

    import sys
    from astropy.table import Table
    from astropy.coordinates import SkyCoord, Angle
    from astropy import units as u
    from regions import CircleSkyRegion
    from get_gsc import get_bbox

    in_file, gsc_file, out_file = sys.argv[1:]
    hdulist, exts, layout = open_with_meta(in_file)
    print "focal plane layout: {0}".format(layout)

    gsc = Table.read(gsc_file, format='ascii.commented_header')
    maglim = 15
    roi = get_bbox(hdulist)
    pad = 60 / 3600.  # 60 arcsec
    ra_pad = pad / np.cos(0.5 * (bs + bn) * np.pi / 180.)
    gsc = gsc[(gsc['ra'] > roi[0] - ra_pad) &
              (gsc['ra'] < roi[1] - ra_pad) &
              (gsc['dec'] > roi[2] - pad) &
              (gsc['dec'] < roi[3] - pad) &
              (gsc['mag'] < maglim)
              ]
    mag = gsc['mag']
    coords = SkyCoord(ra=cat['ALPHA_J2000'], dec=cat['DELTA_J2000'])
    region_string = 'image\n'
    for star in gsc:
        radius = 10 ** ((15 - mag) / 2.5) * 1.4 * u.arcsec
        region_string += 'circle({0}, {1}, {2})'.format(
                star['ra'], star['dec'],)
    for coord in coords:


        region = """
        image
        circle(100, 100, 80)
        box(200, 150, 150, 120, 0)
        """

r = pyregion.parse(region)
mask_1or2 = r.get_mask(shape=(300,300))
        fo.write('fk5;\n')
        for bright_star in gsc:
            bra = bright_star['ra'] * u.degree
            bdec = bright_star['dec'] * u.degree
            bmag = bright_star['mag']
            # if bmag > 90:
            #     bmag = 0.815 * bright_star['JpgMag'] + 1.652
            radius = 10 ** ((15 - bmag) / 2.5) * 1.4 * u.arcsec
            fo.write('circle({0}, {1}, {2}") # color=red\n'.format(
                bra.value, bdec.value, radius.value))
            bcoord = SkyCoord(ra=bra, dec=bdec)
            d2d = bcoord.separation(coord)

    pool = Pool(cpu_count())
    ret = pool.map_async(
            mp_worker,
            [(hdulist[i].data, otaxy) for i, otaxy in exts]
            ).get(9999999)
    for ii, (i, _) in enumerate(exts):
        hdulist[i].data = ret[ii]
    hdulist.writeto(out_file, clobber=True)
