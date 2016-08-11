#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-22 22:52
# Python Version :  2.7.11
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
pp_mosaic.py

Determine and evaluate the zero point
Run swarp to create mosaic
"""


from apus.utils import default_to_main_config
from apus.utils import get_main_config
from apus.utils import tlist_wrapper
from apus.core import dump_config_files
from apus.core import check_config_uptodate
import os


@default_to_main_config
def tlist(inglob, inreg):
    conf = get_main_config()
    # make use of multiple stage of input images
    inglob = ['illum_*.fits', 'masked_*.fits', 'nofri_*.fits']
    t1 = dict(
            name='sex for zeropoint',
            func='sex',
            pipe='transform',
            in_=(inglob, inreg),
            out='sexzp_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.cat',
            dry_run=False,
            params={'CATALOG_TYPE': 'ASCII_HEAD',
                    'DETECT_MINAREA': 3,
                    'DETECT_THRESH': 2,
                    'ANALYSIS_THRESH': 2,
                    'DEBLEND_MINCONT': 0.05,
                    'PHOT_APERTURES': [18, 27, 36, 45, 54, 72, 90, 109],
                    'BACK_SIZE': 128,
                    # 'VERBOSE_TYPE': 'QUIET',
                    },
            outparams=['MAG_APER(8)', 'MAGERR_APER(8)', 'FLUX_MAX',
                       'AWIN_IMAGE', 'BWIN_IMAGE', 'ELONGATION'],
            allow_slice=True,
            )
    t2 = dict(
            name='clean up catalog',
            func='python -u ./recipes/clean_up_catalog.py %s {in} {out}'
            % (conf.band),
            pipe='transform',
            in_=(t1, inreg),
            extras='gsc.cat',  # task 'merge catalog'
            out='cleaned_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.asc',
            dry_run=False,
            follows='merge catalogs',
            )
    t3 = dict(
            name='match to sdss',
            func='./recipes/match_to_sdss.sh {in} {out}',
            pipe='transform',
            in_=(t2, inreg),
            extras='sdss.cat',  # task 'merge catalogs'
            out='matched_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.asc',
            dry_run=False,
            follows='merge catalogs',
            )
    t4 = dict(
            name='create zp catalog',
            func='python -u ./recipes/create_zp_catalog.py {in} {out}',
            pipe='transform',
            in_=(t3, inreg),
            add_inputs='{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.fits',
            out='zp_{imflag[0]}_{id[0]}_{name[0]}_odi_{band[0]}.asc',
            dry_run=False,
            )
    t5 = dict(
            name='create plot list',
            func='python -u ./recipes/create_plot_list.py {in} {out}',
            pipe='collate',
            in_=(t4, inreg),
            out='%s_{band[0]}.plot' % (conf.jobkey),
            dry_run=False,
            )
    t6_prep1 = dict(
            name='create swarp header config',
            func='python -u recipes/create_sc_header.py {out}',
            pipe='originate',
            out=['cat.swarpconf', 'ota.swarpconf', 'none.swarpconf']
            )
    t6_prep2 = dict(
            name='group self calib inputs',
            func='touch {out}',
            pipe='collate',
            in_=(inglob, inreg),
            out='{imflag[0]}_{band[0]}.scgroup',
            follows=conf.tlist[-1],  # from previous TList
            )
    instru = conf.jobkey.split('_')[0]
    t6 = dict(
            name='self calib',
            func='python -u ./recipes/get_flux_scale.py '
                 '{in} {out}',
            pipe='product',
            in_=('*.scgroup',
                 r'(?P<imflag>[^_/]+)_(?P<band>\w)\.scgroup'),
            in2=(['cat.swarpconf', 'ota.swarpconf'],
                 r'(?P<scflag>[^_/]+).swarpconf'),
            replace_inputs=['{imflag[0][0]}_*_odi_{band[0][0]}.fits',
                            'zp_{imflag[0][0]}_*_odi_{band[0][0]}.asc',
                            '{basename[1][0]}{ext[1][0]}'],
            out='%s_{scflag[1][0]}_{imflag[0][0]}_{band[0][0]}.asc' % (instru),
            follows=[t6_prep1, t6_prep2, t4],
            )
    t7 = dict(
            name='plot zp residue',
            func='python -u ./recipes/plot_zp_sc.py {in} {out} save',
            pipe='transform',
            in_="%s_*_?.asc" % (instru, ),
            out='{basename[0]}.eps',
            follows=t6,
            )
    # prepare for product: create dummy file for each swarp input group
    t8_prep = dict(
            name='group swarp inputs',
            func='touch {out}',
            pipe='collate',
            in_=(inglob, inreg),
            out='{imflag[0]}_{name[0]}_odi_{band[0]}.swarpgroup',
            follows=conf.tlist[-1],  # from previous TList
            )
    # three flavors of coadds: sc_per_ota, sc_per_cat, ppa
    # on to different stage of images by imflag
    # conf use relpath due to task default io being jobdir
    confdir = os.path.relpath(conf.confdir, conf.jobdir)
    swarp_conf = os.path.join(confdir, 'conf.swarp_global')
    t8_conf = dict(
            name='create swarp global config',
            func=dump_config_files,
            pipe='originate',
            out=swarp_conf,
            extras=os.path.join(confdir, 'conf.swarp_global.checker'),
            params={
                'PIXELSCALE_TYPE': 'MANUAL',
                'PIXEL_SCALE': 0.20,
                'DELETE_TMPFILES': 'Y',
                'HEADER_SUFFIX': '.none',
                'FSCALE_DEFAULT': 99,
                'BACK_SIZE': 64,
                },
            check_if_uptodate=check_config_uptodate,
            prog='swarp',
            )
    t8 = dict(
            name='swarp',
            func='swarp',
            pipe='product',
            in_=('*.swarpgroup',
                 r'(?P<imflag>[^_/]+)_(?P<name>[^/]+)'
                 r'_odi_(?P<band>\w)\.swarpgroup'),
            in2=('*.swarpconf', r'(?P<scflag>[^_/]+).swarpconf'),
            replace_inputs=['{imflag[0][0]}_*_{name[0][0]}'
                            '_odi_{band[0][0]}.fits',
                            # swarp_conf,
                            # '{basename[1][0]}{ext[1][0]}'
                            ],
            extras=[swarp_conf, '{basename[1][0]}{ext[1][0]}'],
            in_keys=['in', 'conf', 'conf'],
            out=['swarp_{scflag[1][0]}_{imflag[0][0]}_'
                 '{name[0][0]}_odi_{band[0][0]}.fits',
                 'coadd_{scflag[1][0]}_{imflag[0][0]}_'
                 '{name[0][0]}_odi_{band[0][0]}.wht.fits'],
            follows=[t8_prep, t8_conf, t6]
            )
    swarpreg = (r'(?P<prefix>[^_/]+)_'
                r'(?P<scflag>[^_/]+)_'
                r'(?P<imflag>[^_/]+)_(?P<name>[^/]+)_'
                r'odi_(?P<band>\w)\.(?P<ext>.+)')
    t9 = dict(
            name='fix nan pixel',
            func='/home/ma/Codes/apus/scripts/apply_mask.py {in} {out}',
            pipe='transform',
            in_=(t8, swarpreg),
            out='coadd_{scflag[0]}_{imflag[0]}_{name[0]}_odi_{band[0]}.fits'
            )
    # a more general outreg is sufficient for now
    outreg = (r'(?P<prefix>[^_]+_)?coadd_(?P<name>.+)'
              r'_odi_(?P<band>\w)\.fits')
    return tlist_wrapper([t1, t2, t3, t4, t5,
                          t6_prep1, t6_prep2, t6,
                          t7, t8_prep, t8_conf, t8, t9],
                         'coadd_*.fits', outreg)
