#! /bin/sh
#
# check_guider_ota.sh
# Copyright (C) 2016 Jerry Ma <jerry.ma.nk@gmail.com>
#
#----------------------------------------------------------------------------
#"THE BEER-WARE LICENSE" (Revision 42):
#Jerry wrote this file. As long as you retain this notice you
#can do whatever you want with this stuff. If we meet some day, and you think
#this stuff is worth it, you can buy me a beer in return Poul-Henning Kamp
#----------------------------------------------------------------------------



function check() {
    image=$1
    cmd='ds9 -invert -zscale -mosaicimage'
    cmd=$(echo $cmd "${image}")
    echo $cmd
    $cmd
    sleep 1
}

while IFS= read -r entry; do
    image=$(IFS=" " ; set -- $entry ; echo $1)
    if [[ ${image: -5} == ".fits" ]]; then
        check $image
    fi
done
