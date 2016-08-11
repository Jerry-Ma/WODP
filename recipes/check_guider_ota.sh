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
    cmd='ds9 -invert -zscale -tile'
    odi56='10 12 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30'
    for ext in $odi56; do
        cmd=$(echo $cmd "-file" "${image}[${ext}]" "-frame new")
    done
    echo $cmd
    $cmd
}

while IFS= read -r entry; do
    image=$(IFS=" " ; set -- $entry ; echo $1)
    if [[ ${image: -5} == ".fits" ]]; then
        check $image
    fi
done
