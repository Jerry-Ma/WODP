#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-30 20:04
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
pp_common.py
"""


import os
from apus.utils import get_main_config
from apus.utils import default_to_main_config
from apus.utils import tlist_wrapper


def setup_common_config(runkey, band):
    conf = get_main_config()
    conf.jobkey = '_'.join([runkey, band])
    # group by band is probably better than by runkey
    conf.jobdir = os.path.join(band, runkey)
    conf.confdir = os.path.join(conf.jobdir, 'config')
    conf.logdir = 'logs'

    conf.env_overrides = {
            'path_prefix': '/home/ma/Codes/astromatic',
            'scratch_dir': '/run/media/ma/S-Ma/SCRATCH'
            }
    conf.task_io_default_dir = conf.jobdir

    # for masking OTA per image
    conf.ota_flag_file = conf.jobkey + '.ota'
    conf.region_mask_dir = runkey + '_regmask'
    # global naming convention
    conf.inglob = 'ppa_*.fits'
    conf.inreg = (r'(?P<prefix>[^_]+_)?'
                  r'(?P<imflag>[^_]+)_'
                  r'(?P<id>20\d{6}T\d{6}\.\d)_(?P<name>.+)'
                  r'_odi_(?P<band>\w)\.(?P<ext>.+)')
    return conf


@default_to_main_config
def tlist(inglob, inreg):
    conf = get_main_config()
    t1 = dict(
            name='get bright star catalog',
            func='python -u recipes/get_gsc.py {in} {out}',
            pipe='collate',
            in_=(inglob, inreg),
            out='gsc_{name[0]}.cat',
            )
    t2 = dict(
            name='get sdss catalog',
            func='python -u recipes/get_sdss.py {in} {out}',
            pipe='collate',
            in_=(inglob, inreg),
            out='sdss_{name[0]}.cat',
            )
    t3 = dict(
            name='merge catalogs',
            func='python -u recipes/merge_catalogs.py {in} {out}',
            pipe='collate',
            in_=[(t1, r'(?P<kind>sdss|gsc)_.+\.cat'),
                 (t2, r'(?P<kind>sdss|gsc)_.+\.cat')],
            out='{kind[0]}.cat',
            )
    t4 = dict(
            name='apply ota mask',
            func='python -u recipes/apply_ota_mask.py {in} {out}',
            pipe='transform',
            extras=os.path.abspath(conf.ota_flag_file),
            in_=(inglob, inreg),
            out='masked_{id[0]}_{name[0]}_odi_{band[0]}.fits',
            )
    # for quick barging with other tlists
    return tlist_wrapper([t1, t2, t3, t4], 'masked_*.fits', inreg)
