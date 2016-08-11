#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-01 11:23
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
podi_illum.py
"""


from apus.utils import default_to_main_config
from apus.utils import tlist_wrapper


@default_to_main_config
def tlist(inglob, inreg, region_mask_dir, tlist):
    t1 = dict(
            name='sex for illum',
            func='sex',
            pipe='transform',
            in_=(inglob, inreg),
            out_keys=['CATALOG_NAME', 'CHECKIMAGE_NAME'],
            out=['sex_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.cat',
                 'segment_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.fits'],
            params={'CATALOG_TYPE': 'ASCII_HEAD',
                    'DETECT_MINAREA': 3,
                    'DETECT_THRESH': 3,
                    'ANALYSIS_THRESH': 3,
                    'DEBLEND_MINCONT': 0.05,
                    'BACK_SIZE': 128,
                    'CHECKIMAGE_TYPE': 'SEGMENTATION',
                    },
            )
    t2 = dict(
            name='apply object mask',
            func='python -u ./recipes/apply_object_mask.py %s {in} {out}'
            % (region_mask_dir),
            pipe='transform',
            in_=(inglob, inreg),
            add_inputs=t1['out'][1],
            out='sky_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.fits',
            dry_run=False,
            follows=t1,
            )
    t3 = dict(
            name='create illum pattern',
            func='python -u ./recipes/create_illum.py {in} {out}',
            pipe='collate',
            in_=(t2, inreg),
            out='illum_odi_{band[0]}.fits',
            dry_run=False,
            jobs_limit=1,
            )
    t4 = dict(
            name='correct illum',
            func='python -u ./recipes/correct_illum.py {in} {out}',
            pipe='transform',
            in_=(inglob, inreg),
            add_inputs=t3,
            out='illum_{id[0]}_{name[0]}_odi_{band[0]}.fits',
            dry_run=False,
            jobs_limit=1,
            allow_slice=True,
            )
    return tlist_wrapper([t1, t2, t3, t4], 'illum_*.fits', inreg)
