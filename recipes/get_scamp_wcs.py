#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-10 22:25
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
get_scamp_wcs.py
"""


if __name__ == "__main__":
    from astropy.table import Table
    from astropy.io import fits
    import numpy as np
    import sys
    from astropy import wcs

    in_images = sys.argv[1:-2]
    in_file, out_file = sys.argv[-2:]
    tbl = Table.read(in_file, format='ascii.commented_header')
    # get bbox from catalog
    minra = np.min(tbl['ra'])
    maxra = np.max(tbl['ra'])
    mindec = np.min(tbl['dec'])
    maxdec = np.max(tbl['dec'])
    if np.abs(maxra - minra) > 90.:
        raise RuntimeError("catalog spans two wide range {0}--{1}".format(
            maxra, minra))
    cra = (maxra + minra) * 0.5
    cdec = (maxdec + mindec) * 0.5
    width = (maxra - minra) * np.cos(cdec * np.pi / 180.)
    height = (maxdec - mindec)

    hdulist = fits.open(in_images[0])
    ow = wcs.WCS(hdulist[1].header)
    hdulist.close()
    ps = abs(ow.pixel_scale_matrix[0, 0])
    ncol = width / ps
    nrow = height / ps
    ow.wcs.crpix = [ncol / 2, nrow / 2]
    # ow.wcs.cdelt = np.diag(ow.pixel_scale_matrix)
    ow.wcs.crval = [cra, cdec]
    ow.wcs.ctype = ["RA---STG", "DEC--STG"]
    header = ow.to_header()
    # for key in ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2']:
    #     header[key.replace('PC', 'CD')] = header[key]
    header['CD1_1'] = ps
    header['CD2_2'] = -ps
    header['CD1_2'] = 0
    header['CD2_1'] = 0
    keepkeys = ['EQUINOX', 'RADESYS',
                'CTYPE1', 'CUNIT1', 'CRVAL1', 'CRPIX1', 'CD1_1', 'CD1_2',
                'CTYPE2', 'CUNIT2', 'CRVAL2', 'CRPIX2', 'CD2_1', 'CD2_2',
                ]
    for key in header.keys():
        if key not in keepkeys:
            del header[key]
    header.insert(0, ('EXTEND', True))
    header.insert(0, ('NAXIS2', int(nrow + 0.5)))
    header.insert(0, ('NAXIS1', int(ncol + 0.5)))
    header.insert(0, ('NAXIS', 2))
    header.insert(0, ('BITPIX', 0))
    header.insert(0, ('SIMPLE', True))
    print header
    with open(out_file, 'w') as fo:
        fo.write(header.tostring())
