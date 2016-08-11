#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-10 21:36
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
pp_scamp.py
"""


from apus.utils import default_to_main_config
from apus.utils import tlist_wrapper
from apus.utils import get_main_config


@default_to_main_config
def tlist(inglob, inreg):
    conf = get_main_config()
    t1 = dict(
            name='sex for astrometry',
            func='sex',
            pipe='transform',
            in_=(inglob, inreg),
            out='sexwcs_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.fits',
            params=dict(
                CATALOG_TYPE='FITS_LDAC',
                DETECT_MINAREA=5,
                DETECT_THRESH=10.0,
                ),
            )
    t2 = dict(
            name='get scamp refcat wcs',
            func='python -u ./recipes/get_scamp_wcs.py {in} {out}',
            pipe='merge',
            in_=[inglob, 'sdss.cat'],
            out='sdss.wcs',
            follows='merge catalogs'
            )
    t3 = dict(
            name='create scamp refcat',
            func='/home/ma/Codes/apus/scripts/create_refcat.py '
            '{0} {1} {2} {{in}} {{out}}'.format(conf.band, 15, 19),
            pipe='transform',
            in_='sdss.cat',
            out='{basename[0]}.fits',
            follows=t2,
            )
    t4 = dict(
            name='scamp',
            func='scamp',
            pipe='merge',
            in_=[t1, t3],
            in_keys=['in+', 'ASTREFCAT_NAME'],
            params={
                'ASTREF_CATALOG': 'FILE',
                'SOLVE_PHOTOM': 'N',
                'HEADER_SUFFIX': '.hdr_wcs',
                'MOSAIC_TYPE': 'SAME_CRVAL',
                },
            out='mergedref.cat',
            out_keys=['MERGEDOUTCAT_NAME', ],
            )
    return tlist_wrapper([t1, t2, t3, t4],
                         '*.hdr_wcs', inreg)
