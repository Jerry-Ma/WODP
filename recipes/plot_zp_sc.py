#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-25 00:28
# Python Version :  2.7.11
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
plot_zp_sc.py
"""

from __future__ import division
import numpy as np

from pyjerry.mympl import SolarizedColor as sc
from pyjerry.mympl import TexStyle as ts
from pyjerry import mympl

import sys
import re
import os

from astropy.table import Table
import logging
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from mpl_toolkits.axes_grid1 import axes_size
from matplotlib import cm
# from joblib import Parallel
# from joblib import delayed

# import itertools
from mpl_toolkits.axes_grid.inset_locator import inset_axes


def plot_zp_color(ax, bulk, ck):
    ax.set_xlim((-0.2, np.max(bulk['catind']) + 1))
    ax.set_ylim((-0.7, 1.))

    color = bulk[ck['cmag1']] - bulk[ck['cmag2']]
    zp1 = bulk[ck['smag']] - bulk[ck['mag']]
    zp1err = np.hypot(bulk[ck['semag']], bulk[ck['emag']])
    airmass = bulk[ck['airmass']]
    otaxy = bulk[ck['otaxy']]
    catind = bulk['catind']

    pltkw = dict(fmt='o', ms=3, capsize=0, mew=0)
    pltkw_ex = dict(fmt='D', ms=2, capsize=0, mew=1, fillstyle='none')
    # legs = [mympl.get_dummy_leg() for _ in range(len(np.unique(catind)))]
    legs = []

    bota = master['catbota']
    botaunc = master['catbotaunc']
    bcat = master['catbcat']
    bcatunc = master['catbcatunc']
    kins = master['catkins'][0]
    kair = master['catkair'][0]
    bair = master['catbair'][0]
    print "residue: {0} ({1})".format(master['catres'][0],
                                      master['catclipres'][0])
    model = (kins + kair * airmass) * color \
        + bota + bair * airmass + bcat
    res = zp1 - model
    clipmask = master['catclipflag'].astype(dtype=bool)
    # res = res - np.median(res)
    for ii, i in enumerate(np.unique(catind)):
        im = (bulk['catind'] == i) & clipmask
        ex = (bulk['catind'] == i) & (~clipmask)
        # (kins + kair * X) * (color) + bota + bair * X
        if np.any(im):
            leg = ax.errorbar(color[im] + bulk['catind'][im], res[im],
                              yerr=zp1err[im],
                              **pltkw)
            for lc in leg[2]:
                ecolor = lc.get_color()[0]
                eecolor = sc.hsv(ecolor, s=0.2, v=0.9)
                lc.set_color(eecolor)
            shortname = re.match(r'.+_20(?:\d\d)(.+)_odi_\w\..+',
                                 bulk['catfile'][im][0]).group(1).replace(
                                         '_', '-')
            legs.append((leg, ts.tt(shortname)))
            if np.any(ex):
                legex = ax.errorbar(color[ex] + bulk['catind'][ex], res[ex],
                                    yerr=zp1err[ex], color=eecolor, **pltkw_ex)
                for lc in legex[2]:
                    lc.set_color(eecolor)
    # xg = np.linspace(-2, 2, 10)
    ax.axhline(0.05 * np.log10(np.e) * 2.5, ls='--', color=sc.base02)
    ax.axhline(-0.05 * np.log10(np.e) * 2.5, ls='--', color=sc.base02)
    ax.text(-0.1, 0.0, ts.tt(r'\pm 5\%'),
            verticalalignment='center', fontsize=20)

    # inset plot for the best fit parameters
    psize = 200. * mympl.LATEX_INCHES_PER_PT
    px = inset_axes(ax,
                    width=psize,
                    height=psize,
                    loc=1,
                    bbox_to_anchor=(0.05, 0.05, 0.87, 0.9),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                    )
    cpx = inset_axes(ax,
                     width=psize * 0.05,
                     height=psize,
                     loc=1,
                     bbox_to_anchor=(0.05, 0.05, 0.87, 0.9),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
                     )
    # get pval bota array
    bota = np.ones((ck['n_ota_y'], ck['n_ota_x']), dtype='d') * np.nan
    botaunc = np.ones((ck['n_ota_y'], ck['n_ota_x']), dtype='d') * np.nan
    for y in range(*ck['ota_y_range']):
        for x in range(*ck['ota_x_range']):
            if np.any(otaxy == 10 * x + y) > 0:
                bota[y - ck['ota_y_range'][0], x - ck['ota_x_range'][0]] \
                    = master['catbota'][
                                otaxy == 10 * x + y][0]
                botaunc[y - ck['ota_y_range'][0], x - ck['ota_x_range'][0]] \
                    = master['catbotaunc'][
                                otaxy == 10 * x + y][0]
            else:
                print "no data found for OTA {0}{1}".format(x, y)

    pleg = px.imshow(bota, interpolation='none', aspect=1, cmap=cm.coolwarm)
    cb = plt.colorbar(pleg, cax=cpx)
    cb.set_label(ts.tt('b_{OTA}'))
    px.yaxis.set_ticks(range(ck['n_ota_y']))
    px.yaxis.set_ticklabels(['Y{0}'.format(j) for j
                             in reversed(range(*ck['ota_y_range']))])
    px.xaxis.set_ticks(range(ck['n_ota_x']))
    px.xaxis.set_ticklabels(['X{0}'.format(j) for j
                             in range(*ck['ota_x_range'])])
    # plot airmass scatter
    xsize = 200. * mympl.LATEX_INCHES_PER_PT
    xx = inset_axes(ax,
                    width=xsize * 1.5,
                    height=xsize,
                    loc=9,
                    bbox_to_anchor=(0.05, 0.05, 0.87, 0.9),
                    bbox_transform=ax.transAxes,
                    borderpad=0,
                    )
    u_data = []
    for ii, i in enumerate(np.unique(catind)):
        u_x = airmass[bulk['catind'] == i][0]
        u_b = bcat[bulk['catind'] == i][0]
        u_bunc = bcatunc[bulk['catind'] == i][0]
        color = legs[ii][0][0].get_color()
        pltkw_xx = dict(pltkw, ms=9 if ii == 0 else 6, color=color)
        xx.errorbar([u_x], [u_b], yerr=[u_bunc], **pltkw_xx)
        u_data.append((u_x, u_b))
    xx.plot(*zip(*u_data), linestyle='-', color=sc.base1, zorder=0)
    xx.set_xlim((np.min(airmass) - 0.1, np.max(airmass) + 0.1))
    xx.set_xlabel(ts.tt("Airmass"))
    xx.set_ylabel(ts.tt("b_{cat}"))

    # show label text
    label = [
            r'k_{ins}=%.4f' % (kins),
            r'{\langle}b_{ota}{\rangle}=%.4f' % (np.mean(bota)),
            r'{\langle}{\delta}b_{ota}{\rangle}=%.4f' % (np.mean(botaunc)),
            r'{\Delta}b_{ota}=%.4f' % (np.std(bota)),
            r'b_{X}=%.4f' % (bair),
            r'n_{obj}=%d' % (len(zp1)),
            r'{\Delta}res.=%.5f' % (np.std(res)),
            r'n_{obj,clip}=%d' % (len(zp1[clipmask])),
            r'{\Delta}res._{clip}=%.5f' % (np.std(res[clipmask]))
            ]
    ax.text(0.05, 0.92, ts.tt('\n'.join(label)),
            transform=ax.transAxes,
            verticalalignment='top')
    ax.legend(*zip(*legs), loc='lower center', ncol=4)
    return legs


if __name__ == '__main__':

    import get_flux_scale as self_calib
    logger = logging.getLogger(__name__)

    in_file = sys.argv[1]
    if sys.argv[-1] == 'save':
        sys.argv.pop()
        save = True
    else:
        save = False
    plotid, ext = os.path.splitext(os.path.basename(in_file))
    ck = self_calib.get_context(plotid)
    if ext == '.plot':
        cat_list = np.loadtxt(in_file, dtype=str, ndmin=2)
        # create master calib first
        model_flags = sys.argv[2:]
        master = self_calib.get_calibrated_master(
                cat_list[:, 0], ck, model_flags)
        out_file = "fig_{0}.eps".format(plotid)
    elif ext == '.asc':
        # load the master calib
        master = Table.read(in_file, format='ascii.commented_header')
        out_file = sys.argv[2]
    cband = ck['cband']
    # do the plot
    mympl.use_hc_color('kelly')
    canvas = mympl.CanvasOne(
        width=1200,
        aspect=0.6,
        scale=1,
        usetw=False,
        )
    fig, (ax, ) = canvas.parts()

    ax.set_xlabel(ts.tt(r'SDSS {0} - {1} (mag)'.format(*cband)))
    ax.set_ylabel(ts.tt(r'm_{SDSS} - m_{ins} - k_{ins} \times (%s - %s) '
                        r'- b_{OTA} - b_{X}X - b_{cat} (mag)' % cband))

    plot_zp_color(ax, master, ck)
    canvas.save_or_show(out_file,
                        bbox_inches='tight',
                        save=save,
                        )
