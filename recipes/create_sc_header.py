#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Create Date    :  2016-08-04 15:43
# Python Version :  2.7.12
# Git Repo       :  https://github.com/Jerry-Ma
# Email Address  :  jerry.ma.nk@gmail.com
"""
create_sc_header.py
"""


import os


if __name__ == "__main__":
    import sys
    for out_file in sys.argv[1:]:
        header_suffix = os.path.splitext(os.path.basename(out_file))[0]
        with open(out_file, 'w') as fo:
            fo.write("HEADER_SUFFIX .hdr_{0}  #".format(header_suffix))
