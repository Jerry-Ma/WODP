#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-07-31 21:42
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
merge_catalogs.py
"""


from astropy.table import Table
from astropy.table import vstack
from astropy.table import unique


if __name__ == "__main__":
    import sys
    in_files = sys.argv[1:-1]
    out_file = sys.argv[-1]
    tbls = [Table.read(f, format='ascii.commented_header') for f in in_files]
    tbl = unique(vstack(tbls, join_type='exact'))
    tbl.write(out_file, format='ascii.commented_header')
