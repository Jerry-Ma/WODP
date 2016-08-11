#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 14:58
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
create_zp_catalog.py
"""


from astropy.table import Table, Column
from pyjerry.instrument.wiyn import WIYNFact
from common import open_with_meta


if __name__ == "__main__":
    import sys
    from astropy.io import fits
    from astropy.time import Time

    in_file, img_file, out_file = sys.argv[1:]

    renamecol = [
            ('ra', 'SDSS_RA'), ('dec', 'SDSS_DEC'),
            ('u', 'SDSS_MAG_U'), ('err_u', 'SDSS_ERR_U'),
            ('g', 'SDSS_MAG_G'), ('err_g', 'SDSS_ERR_G'),
            ('r', 'SDSS_MAG_R'), ('err_r', 'SDSS_ERR_R'),
            ('i', 'SDSS_MAG_I'), ('err_i', 'SDSS_ERR_I'),
            ('z', 'SDSS_MAG_Z'), ('err_z', 'SDSS_ERR_Z'),
            ('ALPHA_J2000', 'ODI_RA'), ('DELTA_J2000', 'ODI_DEC'),
            ('MAG_AUTO', 'ODI_MAG_AUTO'), ('MAGERR_AUTO', 'ODI_ERR_AUTO'),
            ('XWIN_IMAGE', 'ODI_X'), ('YWIN_IMAGE', 'ODI_Y'),
            ]

    hdulist, exts, layout = open_with_meta(img_file)
    hdulist.close()
    if layout == 'odi56':
        get_ota_xy = WIYNFact.get_ota_xy
    elif layout == 'podi':
        get_ota_xy = WIYNFact.get_ota_xy_podi
    else:
        raise RuntimeError("layout {0} not recognized".format(layout))

    tbl = Table.read(in_file, format='ascii.commented_header')
    for oc, nc in renamecol:
        tbl.rename_column(oc, nc)
    # header column
    hdulist = fits.open(img_file)
    for key in ['AIRMASS', 'EXPMEAS']:
        col = Column([hdulist[0].header[key], ] * len(tbl), name=key)
        tbl.add_column(col)
    # get mjd time
    obstime = Time(hdulist[0].header['DATE-MID'], format='isot', scale='utc')
    col_time = Column([obstime.mjd, ] * len(tbl), name='MJD')
    tbl.add_column(col_time)
    hdulist.close()
    # odixy column
    ota_xy = [get_ota_xy(ext) for ext in tbl['EXT_NUMBER']]
    col_odi = Column(ota_xy, name='ODI_OTA')
    tbl.add_column(col_odi)
    tbl.write(out_file, format='ascii.commented_header')
