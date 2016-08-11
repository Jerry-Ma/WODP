#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 15:46
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
pp_photometry.py

perform photometry and evaluate the depth
"""


from apus.utils import default_to_main_config
from apus.utils import tlist_wrapper


@default_to_main_config
def tlist(inglob, inreg):
    t1 = dict(
            name='do photometry',
            func='sex',
            pipe='transform',
            in_=(inglob, inreg),
            add_inputs='{basename[0]}.wht.fits',
            in_keys=[('in', 'WEIGHT_IMAGE'), ],
            out='sex_{basename[0]}.cat',
            dry_run=False,
            params={'CATALOG_TYPE': 'ASCII_HEAD',
                    'DETECT_MINAREA': 3,
                    'DETECT_THRESH': 1.0,
                    'ANALYSIS_THRESH': 1.0,
                    'DEBLEND_MINCONT': 0.005,
                    'PHOT_APERTURES': [18, 27, 36, 45, 54, 72, 90, 109],
                    'BACK_SIZE': 64,
                    'WEIGHT_TYPE': 'MAP_WEIGHT',
                    # 'WEIGHT_TYPE': 'NONE',
                    'WEIGHT_GAIN': 'Y',
                    'GAIN_KEY': 'GAIN',
                    # 'GAIN': 1.23,
                    'MAG_ZEROPOINT': 25,
                    # 'VERBOSE_TYPE': 'QUIET',
                    },
            outparams=['MAG_APER(8)', 'MAGERR_APER(8)'],
            allow_slice=True,
            )
    t2 = dict(
            name='convert to ascii',
            func='/home/ma/Codes/apus/scripts/sex_to_ascii.py {in} {out}',
            pipe='transform',
            in_=(t1, r'sex_(?P<name>.+).cat'),
            out='phot_{name[0]}.asc',
            )
    t3 = dict(
            name='plot depth histogram',
            func='./recipes/plot_depth_histo.sh 5 {in} {out}',
            pipe='transform',
            in_=t2,
            out='{basename[0]}.png',
            )
    t4 = dict(
            name='create mosaic zp catalog',
            func='./recipes/match_to_sdss.sh {in} {out}',
            pipe='transform',
            in_=(t2, 'phot_(?P<name>.+)\.asc'),
            extras='sdss.cat',
            out='zp_{name[0]}.asc',
            follows='merge catalogs',
            )
    t5 = dict(
            name='plot mosaic zp',
            func='python -u recipes/plot_phot_zp.py {in} {out} 1',
            pipe='transform',
            in_=t4,
            out='{basename[0]}.eps',
            )

    return tlist_wrapper([
        t1, t3, t4, t5, t2], 'phot_*.asc', 'phot_(?P<name>.+).asc')
