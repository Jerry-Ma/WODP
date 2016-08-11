#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-24 18:05
# Python Version :  2.7.11
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
get_sdss.py
"""

from get_gsc import get_bbox
from get_gsc import to_ds9_box
from get_gsc import merge_bbox


if __name__ == "__main__":
    import os
    import re
    import sys
    from astroquery.sdss import SDSS

    in_files = sys.argv[2:-1]
    out_file = sys.argv[-1]
    out_reg = os.path.splitext(out_file)[0] + '.reg'
    sdss_filter = re.match(r'.+odi_(\w).+', in_files[0]).group(1)
    box = None
    with open(out_reg, 'w') as fo:
        fo.write('global color=red\n')
        for image in in_files:
            print "working on image:", image
            ibox = get_bbox(image)
            icra, icdec, iwidth, iheight = to_ds9_box(ibox)
            fo.write('fk5; box({0},{1},{2},{3}, 0)\n'.format(
                    icra, icdec, iwidth, iheight))
            box = merge_bbox(box, ibox)
        cra, cdec, width, height = to_ds9_box(box)
        fo.write('fk5; box({0},{1},{2},{3}, 0)'.format(
                cra, cdec, width, height))

    # query sdss and get reference catalog
    sql_query = [
        "SELECT ra,dec,raErr,decErr,u,err_u,g,err_g,r,err_r,i,err_i,z,err_z",
        "FROM Star WHERE",
        "    ra BETWEEN {min_ra:f} and {max_ra:f}",
        "AND dec BETWEEN {min_dec:f} and {max_dec:f}",
        "AND ((flags_{filter:s} & 0x10000000) != 0)",     # detected in BINNED1
        # not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP,
        # SATURATED, or BAD_COUNTS_ERROR"
        "AND ((flags_{filter:s} & 0x8100000c00a4) = 0)",
        # not DEBLEND_NOPEAK or small PSF error
        "AND (((flags_{filter:s} & 0x400000000000) = 0) or "
        "(psfmagerr_{filter:s} <= 0.2))",
        # not INTERP_CENTER or not COSMIC_RAY
        "AND (((flags_{filter:s} & 0x100000000000) = 0) or "
        "(flags_{filter:s} & 0x1000) = 0)"]
    sql_query = '\n'.join(sql_query).format(
            filter=sdss_filter,
            min_ra=box[0], max_ra=box[1],
            min_dec=box[2], max_dec=box[3]
            )
    print sql_query
    stdstar = SDSS.query_sql(sql_query)
    stdstar.write(out_file, format='ascii.commented_header')
