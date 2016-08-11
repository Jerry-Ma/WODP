#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-16 20:14
# Python Version :  2.7.11
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
get_gsc.py

Get GSC bright star catalog for bright star rejection purpose
"""


from pyjerry.instrument.wiyn import WIYNFact
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np


def get_bbox(image):
    hdulist = fits.open(image)
    coord = SkyCoord(hdulist[0].header['RA'], hdulist[0].header['DEC'],
                     unit=(u.hourangle, u.degree))
    hdulist.close()
    ra = coord.ra.degree
    dec = coord.dec.degree
    (w, e), (s, n) = WIYNFact.get_bbox()
    bw = ra + w / np.cos(dec * np.pi / 180.)
    be = ra + e / np.cos(dec * np.pi / 180.)
    bs = dec + s
    bn = dec + n
    return bw, be, bs, bn


def get_min_max(l):
    minl = np.inf
    maxl = -np.inf
    dist = 0
    for i in l:
        for j in l:
            # TODO get handle of degree wrapping
            if np.abs(i - j) > dist:
                minl = min(i, j)
                maxl = max(i, j)
                dist = np.abs(i - j)
    if dist > 180.:
        print l
        raise RuntimeError(
            "there is degree wrapping that that code is not able to handle")
    return minl, maxl


def merge_bbox(box1, box2):
    if box1 is None:
        return box2
    elif box2 is None:
        return box1
    else:
        l1, r1, b1, t1 = box1
        l2, r2, b2, t2 = box2
        ramin, ramax = get_min_max([l1, l2, r1, r2])
        decmin, decmax = get_min_max([b1, b2, t1, t2])
        return ramin, ramax, decmin, decmax


def to_ds9_box(box):
    cra = (box[0] + box[1]) * 0.5
    cdec = (box[2] + box[3]) * 0.5
    dra = (box[1] - box[0]) * np.cos(cdec * np.pi / 180.)
    ddec = (box[3] - box[2])
    return cra, cdec, dra, ddec


if __name__ == '__main__':
    import os
    import sys
    in_files = sys.argv[1:-1]
    out_file = sys.argv[-1]
    out_reg = os.path.splitext(out_file)[0] + '.reg'

    box = None
    with open(out_reg, 'w') as fo:
        fo.write('global color=red\n')
        for image in in_files:
            print "working on image:", image
            ibox = get_bbox(image)
            icra, icdec, iwidth, iheight = to_ds9_box(ibox)
            fo.write('fk5; box({0},{1},{2},{3}, 0)\n'.format(
                    icra, icdec, iwidth, iheight))
            box = merge_bbox(box, ibox)
        cra, cdec, width, height = to_ds9_box(box)
        fo.write('fk5; box({0},{1},{2},{3}, 0)'.format(
                cra, cdec, width, height))

    # query guide star catalog
    url = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx'
    payload = {
            'BBOX': '{0},{1},{2},{3}'.format(box[0], box[2], box[1], box[3]),
            'FORMAT': 'VOTable', 'CAT': 'GSC23'}
    import requests
    from astropy.table import Table
    from astropy.table import Column
    response = requests.get(url=url, params=payload)
    if response.status_code == requests.codes.ok:
        outcol = ['ra', 'dec', 'FpgMag', 'JpgMag', 'NpgMag']
        with open(out_file, 'w') as fo:
            fo.write(response.text)
        # filter the catalog
        cat = Table.read(out_file, format='votable')
        cat = cat[(cat['FpgMag'] < 16) |
                  (cat['JpgMag'] < 16) |
                  (cat['NpgMag'] < 16)
                  ][outcol]
        # create synthesis mag column based on mean color
        mcol = 'FpgMag'
        mag = cat[mcol]
        for col in ['JpgMag', 'NpgMag']:
            gmask = (mag < 90) & (cat[col] < 90)
            offset = np.mean(mag[gmask] - cat[col][gmask])
            bmask = (mag > 90) & (cat[col] < 90)
            mag[bmask] = cat[col][bmask] + offset
        cat.add_column(Column(mag, name='mag'))
        cat = cat[cat['mag'] < 15]
        cat.write(out_file, format='ascii.commented_header')
