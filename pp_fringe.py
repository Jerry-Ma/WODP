#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-30 21:34
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
pp_fringe.py
"""


from apus.utils import default_to_main_config
from apus.utils import tlist_wrapper


@default_to_main_config
def tlist(inglob, inreg):
    t1 = dict(
            name='create fringe template',
            func='python -u ./recipes/create_fringe.py {in} {out}',
            pipe='collate',
            in_=(inglob, inreg),
            out='fringe_odi_{band[0]}.fits',
            dry_run=False,
            jobs_limit=1,
            )
    t2 = dict(
            name='decompose fringe template',
            func='python -u ./recipes/decompose_fringe.py {in} {out}',
            pipe='transform',
            in_=(t1, r'(?P<prefix>[^_]+)_odi_(?P<band>\w)\.fits'),
            out=['hf_odi_{band[0]}.fits', 'lf_odi_{band[0]}.fits'],
            dry_run=False,
            jobs_limit=1,
            )
    t3 = dict(
            name='remove fringe',
            func='python -u ./recipes/remove_fringe.py {in} {out}',
            pipe='transform',
            in_=(inglob, inreg),
            add_inputs=[t1, t2['out'][0]],
            out='nofri_{id[0]}_{name[0]}_odi_{band[0]}.fits',
            dry_run=False,
            jobs_limit=1,
            allow_slice=True,
            follows=t2,
            )
    return tlist_wrapper([t1, t2, t3], 'nofri_*.fits', inreg)
