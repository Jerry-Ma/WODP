#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 16:38
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
get_flux_scale.py

From a list of image + zp catalog, perform self calibration
and generate swarp header files
"""


from __future__ import division
import numpy as np
from pyjerry.instrument.wiyn import WIYNFact

import sys
import os

from scipy.stats import sigmaclip
from astropy.table import Table
from astropy.table import vstack
from astropy.table import Column

import logging
import lmfit
# from joblib import Parallel
# from joblib import delayed

# import itertools
import time


def get_master_calib(tablelist, ck):
    '''get calibration table files from input file list'''
    logger = logging.getLogger(__name__)
    # load table and compile to a master calib table
    logger.info('create master')
    master = []
    for i, fname in enumerate(tablelist):
        try:
            tbl = Table.read(fname, format='ascii.commented_header')
        except IOError:
            logger.warning("{0} does not exist".format(fname))
            continue
        # reject low exposure time
        exptime = np.mean(tbl['EXPMEAS'])
        if exptime < 100.:
            logger.warning("{0} discarded due to short exptime {1}"
                           .format(fname, exptime))
            continue
        # global zero-th order zp by a 5 sigma clip
        zp = tbl[ck['smag']] - tbl[ck['mag']]
        clipped, lolim, uplim = sigmaclip(zp, ck['clip0'], ck['clip0'])
        zp0 = np.median(clipped)
        logger.info('zero-th order zp clipping ({0:.1f} sigma):'
                    ' {1:+.3f} {2:.3f} {3:+.3f}'
                    .format(ck['clip0'], lolim - zp0, zp0, uplim - zp0))
        smag_lo, smag_up = ck['smaglims']
        qualitymask = (smag_lo < tbl[ck['smag']]) & \
                      (tbl[ck['smag']] < smag_up) & \
                      (zp > lolim) & (zp < uplim)
        tbl = tbl[qualitymask]
        n = len(tbl)
        if n < 10:
            logger.warning("{0} contains too few objects".format(fname))
            continue
        logger.info('read in {0} ({1})'.format(fname, n))
        # append the table level info
        col_ind = Column([i] * n, name='catind')
        col_cfn = Column([os.path.basename(fname)] * n, name='catfile')
        col_zp0 = Column([zp0] * n, name='catzp0')
        for col in [col_ind, col_cfn, col_zp0]:
            tbl.add_column(col)
        # place holder for the self calibration result
        # (kins + kair * X) * (color) + bota + bair * X + bcat
        for name in [
                'catkins', 'catkinsunc', 'catkair', 'catkairunc',
                'catbota', 'catbotaunc', 'catbair', 'catbairunc',
                'catbcat', 'catbcatunc', 'catclipflag',
                'catres', 'catclipres'
                ]:
            col = Column([0.0] * n, name=name)
            tbl.add_column(col)
        # col_kins = Column([0.0] * n, name='catkins')
        # col_kair = Column([0.0] * n, name='catkair')
        # col_bota = Column([0.0] * n, name='catbota')
        # col_bair = Column([0.0] * n, name='catbair')
        # col_bcat = Column([0.0] * n, name='catbcat')
        # for col in [col_ind, col_cfn, col_zp0,
        #             col_kins, col_kair, col_bota, col_bair, col_bcat]:
        #     tbl.add_column(col)
        master.append(tbl)
    master = vstack(master, join_type='exact')
    return master


def sc_func(params, color, zp, zperr, otaxy, catind, airmass):
    # (kins + kair * X) * (color) + bota + bair * X + bcat
    bota = np.array([params['b{0:.0f}'.format(o)].value for o in otaxy])
    bcat = np.array([params['bcat{0:.0f}'.format(o)].value for o in catind])
    model = (params['kins'].value + params['kair'].value * airmass) * color \
        + bota + params['bair'].value * airmass + bcat
    return (model - zp) / zperr


def self_calibrate(bulk, ck, model_flags=('color', 'ota', 'cat', 'airmass')):
    '''
    For each exposure/ota, fit an offset so that the over all
    dispersion is minimized
    '''
    color = bulk[ck['cmag1']] - bulk[ck['cmag2']]
    zp1 = bulk[ck['smag']] - bulk[ck['mag']]
    zp1err = np.hypot(bulk[ck['semag']], bulk[ck['emag']])
    params = lmfit.Parameters()
    params.add('kins', value=ck['kins'], vary='color' in model_flags)
    params.add('kair', value=0., vary=False)
    params.add('bair', value=0, vary='airmass' in model_flags)
    otaxy = bulk[ck['otaxy']]
    catind = bulk['catind']
    airmass = bulk[ck['airmass']]
    for xy in np.unique(otaxy):
        vary = False if xy == 33 else 'ota' in model_flags
        params.add('b{0:.0f}'.format(xy), value=np.mean(zp1), vary=vary)
    for ind in np.unique(catind):
        vary = 'cat' in model_flags
        params.add('bcat{0:.0f}'.format(ind), 0.0, vary=vary)
    outparams = lmfit.minimize(
            sc_func, params, args=(color, zp1, zp1err,
                                   otaxy, catind, airmass)
            ).params
    # do a sigma clipping on the de-trended data
    residue = sc_func(
            outparams, color, zp1, zp1err, otaxy, catind, airmass) * zp1err
    clipped, lo, up = sigmaclip(residue, ck['clipsc'], ck['clipsc'])
    clipmask = (residue >= lo) & (residue <= up)
    # do the fitting again
    outparams = lmfit.minimize(
            sc_func, outparams,
            args=(color[clipmask], zp1[clipmask], zp1err[clipmask],
                  otaxy[clipmask], catind[clipmask],
                  airmass[clipmask])).params

    parvals = outparams.valuesdict()
    paruncs = {k: outparams[k].stderr for k in parvals.keys()}
    return parvals, paruncs, clipmask, np.std(residue), np.std(clipped)


def get_calibrated_master(in_cats, ck, model_flags):
    logger = logging.getLogger(__name__)
    logger.info("run self calibrate with plotid {0}".format(ck['plotid']))
    master = get_master_calib(in_cats, ck)
    logger.info("self calibration with flags: {0}".format(model_flags))
    start_time = time.time()
    sc_ret = self_calibrate(master, ck, model_flags=model_flags)
    logger.info("self calibration finished after {0}s"
                .format(time.time() - start_time))
    # save the data to master_calib
    outpvals, outpuncs, clipmask, rstd, rcstd = sc_ret
    master['catkins'] = outpvals['kins']
    master['catkinsunc'] = outpuncs['kins']
    master['catkair'] = outpvals['kair']
    master['catkairunc'] = outpuncs['kair']
    master['catbota'] = np.array([outpvals['b{0:.0f}'.format(o)]
                                  for o in master[ck['otaxy']]])
    master['catbotaunc'] = np.array([outpuncs['b{0:.0f}'.format(o)]
                                     for o in master[ck['otaxy']]])
    master['catbair'] = outpvals['bair']
    master['catbairunc'] = outpuncs['bair']
    master['catbcat'] = np.array([outpvals['bcat{0:.0f}'.format(o)]
                                  for o in master['catind']])
    master['catbcatunc'] = np.array([outpuncs['bcat{0:.0f}'.format(o)]
                                     for o in master['catind']])
    master['catclipflag'][clipmask] = 1
    master['catres'] = rstd
    master['catclipres'] = rcstd
    return master


def get_context(plotid):

    jobkey, band = plotid.rsplit('_', 1)
    if jobkey.startswith('podi'):
        ota_x_range = [2, 5]
        ota_y_range = [2, 5]
        layout = 'podi'
    elif jobkey.startswith('odi'):
        ota_x_range = [1, 6]
        ota_y_range = [1, 7]
        layout = 'odi56'

    n_ota_x = ota_x_range[1] - ota_x_range[0]
    n_ota_y = ota_y_range[1] - ota_y_range[0]

    cband = {'u': ('u', 'g'),
             'g': ('g', 'r'),
             'r': ('r', 'i'),
             'i': ('r', 'i'),
             'z': ('i', 'z'),
             }
    ck = {
        'smag': 'SDSS_MAG_{0}'.format(band.upper()),
        'semag': 'SDSS_ERR_{0}'.format(band.upper()),
        'smaglims': dict(z=(0, 19), u=(0, 19)).get(band, (0, 21)),
        # 'mag': 'ODI_MAG_AUTO',
        # 'emag': 'ODI_ERR_AUTO',
        'mag': 'MAG_APER_3',
        'emag': 'MAGERR_APER_3',
        'otaxy': 'ODI_OTA',
        'plotid': plotid,
        'jobkey': jobkey,
        'layout': layout,
        'band': band,
        'cband': cband[band],
        'cmag1': 'SDSS_MAG_{0}'.format(cband[band][0].upper()),
        'cmag2': 'SDSS_MAG_{0}'.format(cband[band][1].upper()),
        'clip0': 5.0,
        'clipsc': 3.0,
        'airmass': 'AIRMASS',
        'kins': dict(g=0.14, z=-0.1342).get(band, 0.),
        'n_ota_x': n_ota_x,
        'n_ota_y': n_ota_y,
        'ota_x_range': ota_x_range,
        'ota_y_range': ota_y_range,
        }
    return ck


def look_up_images(images, cat):
    bname = os.path.splitext(os.path.basename(cat))[0].split('_', 1)[-1]
    _candidate = [f for f in images if bname in f]
    if len(_candidate) == 1:
        return _candidate[0]
    else:
        raise RuntimeError("unable to find image for catalog {0}".format(cat))


if __name__ == "__main__":
    logging.basicConfig(format='[%(name)s] %(message)s', level=logging.DEBUG)
    # model_flags = sys.argv[1].split(',')
    # header_suffix = sys.argv[2]
    # in_images = sys.argv[3:-1:2]
    # in_cats = sys.argv[4:-1:2]
    # out_master = sys.argv[-1]
    in_files = sys.argv[1:-2]
    insize = int(len(in_files) / 2)
    in_images = in_files[:insize]
    in_cats = in_files[insize:]
    in_header = sys.argv[-2]
    out_master = sys.argv[-1]
    model_flags_key = os.path.splitext(os.path.basename(in_header))[0]
    model_flags = dict(ota='color,ota,cat', cat='color,cat')[model_flags_key]
    with open(in_header, 'r') as fo:
        for ln in fo.readlines():
            if ln.startswith("HEADER_SUFFIX"):
                header_suffix = ln.split()[1]
                break
        else:
            raise RuntimeError("no keyword found for header suffix")
    plotid = os.path.splitext(os.path.basename(out_master))[0]
    ck = get_context(plotid)
    master = get_calibrated_master(in_cats, ck, model_flags)
    if ck['layout'] == 'odi56':
        get_ota_xy = WIYNFact.get_ota_xy
    elif ck['layout'] == 'podi':
        get_ota_xy = WIYNFact.get_ota_xy_podi
    else:
        raise RuntimeError('layout {0} not recognized'.format(ck['layout']))
    # generate fluxscale header
    catind = np.unique(master['catind'])
    incat_basenames = map(os.path.basename, in_cats)
    for i in catind:
        subbulk = master[master['catind'] == i]
        catfile = subbulk['catfile'][0]
        airmass = subbulk[ck['airmass']][0]
        # figure out hdrfile name from catfile through indexing
        hdrfile = os.path.join(
            os.path.dirname(out_master),
            os.path.basename(look_up_images(in_images, catfile)
                             ).replace('.fits', header_suffix)
            )
        print "processing swarp header {0}".format(hdrfile)
        with open(hdrfile, 'w') as fo:
            for ext in range(1, ck['n_ota_x'] * ck['n_ota_y'] + 1):
                ota_xy = get_ota_xy(ext)
                extentry = subbulk[subbulk[ck['otaxy']] == ota_xy]
                if len(extentry) == 0:
                    k = b = s = 99
                else:
                    extentry = extentry[0]
                    k = extentry['catkins']
                    b = extentry['catbota'] + extentry['catbair'] * airmass + \
                        extentry['catbcat']
                    s = str(10 ** ((25.0 - b) / 2.5))
                fo.write("{0:8s}={1:>22s}{2:s}\n".format(
                    "MYCOLOR", str(k), " / color term {0}-{1}"
                    .format(*ck['cband'])))
                fo.write("{0:8s}={1:>22s}{2:s}\n".format(
                    "MYZEROP", str(b), " / self-calibrated zero point"))
                fo.write("{0:8s}={1:>22s}{2:s}\n".format(
                    "FLXSCALE", str(s), " / sc flux scale to zp=25"))
                fo.write("END     \n")
    master.write(out_master, format='ascii.commented_header')
    print "master calib file saved: {0}".format(out_master)
