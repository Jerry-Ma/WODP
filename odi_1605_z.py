#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-30 19:35
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
odi_1605_z.py

Data obtained during May, 2016.
The filter used is only z band
There are fringes from the PPA reduced images

The script is an Apus job factory, i.e.,
pipeline will be created on the fly with parameters
such as band, root dir, jobkey, etc.

Steps:
    0. prepare
    1. fringe
    2. zeropoint
    3. photometry

"""

import os
import pp_common
import pp_fringe
# import pp_illum
import pp_scamp
import pp_mosaic
import pp_photometry
from apus import core

runkey, band = 'odi_1605', 'z'
input_root = '/run/media/ma/S-Ma/WIYN/May16_ppa/calibrated'
input_glob = '20??????T*/*odi_{0}.????.fits'.format(band)

conf = pp_common.setup_common_config(runkey, band)
conf.inputs = (
    os.path.join(input_root, input_glob),
    (r'.+(?P<id>20\d{6}T\d{6}\.\d)_(?P<name>\w+)_odi_(?P<band>\w)'
     r'\.(?P<qr>\d{4})/.+\.fits'),
    os.path.join(conf.jobdir, 'ppa_{id[0]}_{name[0]}_odi_{band[0]}.fits')
    )
pp_common.tlist()               \
    .chain(pp_fringe.tlist)     \
    .chain(pp_scamp.tlist)      \
    .chain(pp_mosaic.tlist)     \
    .chain(pp_photometry.tlist) \

core.bootstrap()
