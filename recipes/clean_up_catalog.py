#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 14:41
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
clean_up_catalog.py

read the bright star catalog, calculate mask radius, and
remove objects that falls within it
also remove edge sources
"""


import os
import re
import sys
import itertools
from pyjerry.instrument import wiyn


def get_edgemask(xs, ys, e):
    wl = wiyn.WIYNLayout(binning=1)
    mask = np.zeros_like(xs, dtype=bool)
    for cj, ci in itertools.product(range(wl.NCX), range(wl.NCY)):
        (cl, cr), (cb, ct) = wl.get_cell_rect(0, 0, cj, ci)
        mask = mask | ((xs > cl + e) & (xs < cr - e) &
                       (ys > cb + e) & (ys < ct - e))
    return mask


def get_bpmask(xs, ys, exts, pad=30):
    if np.max(exts) > 15:  # odi56
        layout = 'odi56'
        get_ota_xy = wiyn.WIYNFact.get_ota_xy
        ota_order = wiyn.WIYNFact.ota_order
    else:
        layout = 'podi'
        get_ota_xy = wiyn.WIYNFact.get_ota_xy_podi
        ota_order = wiyn.WIYNFact.ota_order_podi
    print "focal plane layout {0}".format(layout)
    if layout == 'odi56':
        bpmdir = '/home/ma/Codes/podi/trunk/.bpm/odi_5x6'
    elif layout == 'podi':
        bpmdir = '/home/ma/Codes/podi/trunk/.bpm/podi/'
    else:
        raise RuntimeError('focal plane layout {0} not recognized'
                           .format(layout))
    mask = np.ones_like(xs, dtype=bool)
    bpms = {}
    for otaxy in ota_order:
        bpm = []
        bpm_file = os.path.join(bpmdir, 'bpm_xy{0}.reg'.format(otaxy))
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
        bpms[otaxy] = bpm
    for i in range(len(xs)):
        otaxy = get_ota_xy(exts[i])
        bpm = bpms[otaxy]
        for l, r, b, t in bpm:
            if xs[i] > l - pad and xs[i] < r + pad and \
               ys[i] > b - pad and ys[i] < t + pad:
                mask[i] = False
    return mask


if __name__ == "__main__":

    import numpy as np
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.table import Table

    band = sys.argv[1]
    mmag = dict(g=15, z=15, u=12).get(band)
    in_file = sys.argv[2]
    gsc_file = sys.argv[3]
    out_file = sys.argv[4]
    out_reg = os.path.splitext(out_file)[0] + '.reg'

    pad = 20  # number of pixels to edge
    cat = Table.read(in_file, format='ascii.sextractor')
    # remove edge sources and bpm sources
    edgemask = get_edgemask(
            cat['XWIN_IMAGE'], cat['YWIN_IMAGE'], pad)
    bpmask = get_bpmask(
            cat['XWIN_IMAGE'], cat['YWIN_IMAGE'], cat['EXT_NUMBER'], pad=pad)
    cat = cat[edgemask & bpmask]

    # handle bright sources
    gsc = Table.read(gsc_file, format='ascii.commented_header')
    # trim gsc first
    roi = (np.min(cat['ALPHA_J2000']), np.max(cat['ALPHA_J2000']),
           np.min(cat['DELTA_J2000']), np.max(cat['DELTA_J2000']))
    roipad = 20 / 3600.  # 20 arcsec
    gsc = gsc[(gsc['ra'] > roi[0] - pad / 15.) &
              (gsc['ra'] < roi[1] + pad / 15.) &
              (gsc['dec'] > roi[2] - pad) &
              (gsc['dec'] < roi[3] + pad)
              ]
    # remove bright sources
    coord = SkyCoord(ra=cat['ALPHA_J2000'], dec=cat['DELTA_J2000'])
    # loop over gsc and get matches and remove
    goodmask = np.ones((len(cat), ), dtype=bool)
    with open(out_reg, 'w') as fo:
        fo.write('fk5;\n')
        for bright_star in gsc:
            bra = bright_star['ra'] * u.degree
            bdec = bright_star['dec'] * u.degree
            bmag = bright_star['mag']
            # if bmag > 90:
            #     bmag = 0.815 * bright_star['JpgMag'] + 1.652
            radius = 10 ** ((mmag - bmag) / 2.5) * 1.4 * u.arcsec
            fo.write('circle({0}, {1}, {2}") # color=red\n'.format(
                bra.value, bdec.value, radius.value))
            bcoord = SkyCoord(ra=bra, dec=bdec)
            d2d = bcoord.separation(coord)
            goodmask = goodmask & (d2d > radius)
    cat = cat[goodmask]
    cat.write(out_file, format='ascii.commented_header')
