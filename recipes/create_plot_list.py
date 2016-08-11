#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 15:12
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
create_plot_list.py
"""


if __name__ == "__main__":
    import sys
    from astropy.table import Table

    in_files = sys.argv[1:-1]
    out_file = sys.argv[-1]

    data = []
    for fname in in_files:
        tbl = Table.read(fname, format='ascii.commented_header')
        mjd = tbl['MJD'][0]
        data.append((fname, mjd))
    data.sort(key=lambda x: x[1])
    with open(out_file, 'w') as fo:
        fo.write("# catalog mjd\n")
        for d in data:
            fo.write("{0} {1}\n".format(*d))
