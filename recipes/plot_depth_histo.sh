#! /bin/sh
#
# plot_depth_histo.sh
# Copyright (C) 2016 Jerry Ma <jerry.ma.nk@gmail.com>
#
#----------------------------------------------------------------------------
#"THE BEER-WARE LICENSE" (Revision 42):
#Jerry wrote this file. As long as you retain this notice you
#can do whatever you want with this stuff. If we meet some day, and you think
#this stuff is worth it, you can buy me a beer in return Poul-Henning Kamp
#----------------------------------------------------------------------------

snr=$1
cat=$2
out=$3

title=${out%.*}

stilts plot2plane \
    xmin=20 xmax=30 xlabel='MAG_AUTO (mag)' ylabel='Count' binsize=0.1 \
    legpos=0.95,0.95 title=$title \
    in=$2 ifmt=ascii x=mag_auto\
    layer1=histogram color1=red leglabel1='All' \
    layer2=histogram color2=blue leglabel2="SNR>${snr}" \
    icmd2="select \"MAGERR_AUTO<1.0857/${snr}\"" \
    out=$3
