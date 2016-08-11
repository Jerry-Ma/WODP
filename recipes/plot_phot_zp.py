#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-05 12:19
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
plot_phot_zp.py
"""

from __future__ import division
import numpy as np

from pyjerry.mympl import SolarizedColor as sc
from pyjerry.mympl import TexStyle as ts
from pyjerry import mympl
from scipy.stats import sigmaclip


def plot_zp_color(ax, bulk, ck):
    bulk = bulk[(bulk[ck['smag']] > ck['smaglims'][0]) &
                (bulk[ck['smag']] < ck['smaglims'][1])]
    color = bulk[ck['cmag1']] - bulk[ck['cmag2']]
    zp = bulk[ck['smag']] - bulk[ck['mag']] - ck['zp'] - ck['kins'] * color
    zp0 = np.median(zp)
    print "median residue", np.median(zp)
    if ck['scflag'] == 'oditool':
        print "adjust zp level for odi-tools"
        zp = zp - zp0
    zperr = np.hypot(bulk[ck['semag']], bulk[ck['emag']])
    clipped, lo, up = sigmaclip(zp, ck['clip'], ck['clip'])
    im = (zp < up) & (zp > lo)
    ex = ~im
    ax.set_xlim((-0.5, 2))
    ax.set_ylim((-0.7, 1.))

    pltkw = dict(ms=6, fmt='o')
    pltkw = dict(fmt='o', ms=3, capsize=0, mew=0)
    pltkw_ex = dict(fmt='D', ms=2, capsize=0, mew=1, fillstyle='none')
    legs = []
    if np.any(im):
        leg = ax.errorbar(color[im], zp[im], yerr=zperr[im], **pltkw)
        for lc in leg[2]:
            ecolor = lc.get_color()[0]
            eecolor = sc.hsv(ecolor, s=0.2, v=0.9)
            lc.set_color(eecolor)
        legs.append((leg, ts.tt('clipped')))
        if np.any(ex):
            legex = ax.errorbar(color[ex], zp[ex], yerr=zperr[ex],
                                color=eecolor, **pltkw_ex)
            for lc in legex[2]:
                lc.set_color(eecolor)
    ax.axhline(0.05 * np.log10(np.e) * 2.5, ls='--', color=sc.base02)
    ax.axhline(-0.05 * np.log10(np.e) * 2.5, ls='--', color=sc.base02)
    ax.text(-0.1, 0.0, ts.tt(r'\pm 5\%'),
            verticalalignment='center', fontsize=20)
    label = [
            r'{\Delta}res.=%.5f' % (np.std(zp)),
            r'n_{obj,clip}=%d' % (len(zp[im])),
            r'{\Delta}res._{clip}=%.5f' % (np.std(zp[im]))
            ]
    ax.text(0.05, 0.92, ts.tt('\n'.join(label)),
            transform=ax.transAxes,
            verticalalignment='top')


if __name__ == '__main__':

    import sys
    from astropy.table import Table
    import re
    import glob
    import os

    cat_file, out_name = sys.argv[1:3]
    cat = Table.read(cat_file, format='ascii.commented_header')
    match = re.match('([^/]+)_([ugriz]).asc', os.path.basename(cat_file))
    plotid = match.group(1)
    band = match.group(2)
    # figure out color term from hdr file
    print plotid
    match = re.match(r'zp_coadd_(?P<scflag>[a-z]+)_(?P<imflag>[a-z]+)'
                     r'_(?P<name>.+)_odi',
                     plotid)
    scflag = match.group('scflag')
    imflag = match.group('imflag')
    objname = match.group('name')
    if scflag in ['none', 'oditool']:
        kins = dict(g=0.14, ).get(band, 0.)
        print "default color term: {0}".format(kins)
    else:
        hdr_file = glob.glob(
            os.path.join(
                os.path.dirname(cat_file),
                '{0}_*_{1}_odi_{2}.hdr_{3}'.format(imflag, objname,
                                                   band, scflag)
                                        ))[0]
        with open(hdr_file, 'r') as fo:
            for ln in fo.readlines():
                if ln.strip().startswith('MYCOLOR'):
                    kins = float(ln.split()[2].strip())
                    print "color term from hdr: {0}".format(kins)
    cband = {'u': ('u', 'g'),
             'g': ('g', 'r'),
             'r': ('r', 'i'),
             'i': ('r', 'i'),
             'z': ('i', 'z'),
             }
    ck = {
        'smag': band,
        'semag': 'err_{0}'.format(band),
        'smaglims': dict(z=(0, 19), ).get(band, (0, 21)),
        'mag': 'MAG_AUTO',
        'emag': 'MAGERR_AUTO',
        # 'mag': 'MAG_APER_3',
        # 'emag': 'MAGERR_APER_3',
        'plotid': plotid,
        'band': band,
        'cband': cband[band],
        'cmag1': cband[band][0],
        'cmag2': cband[band][1],
        'clip': 3.0,
        'zp': 0.,
        'kins': kins,
        'scflag': scflag
        }

    canvas = mympl.CanvasOne(
        width=800,
        aspect=0.618,
        scale=1,
        usetw=False,
        )
    fig, (ax, ) = canvas.parts()
    plot_zp_color(ax, cat, ck)

    canvas.save_or_show(out_name,
                        bbox_inches='tight',
                        # pad_inches=0,
                        )
