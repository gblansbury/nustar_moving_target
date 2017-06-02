#!/bin/bash

basedir=/home/nustar1/nustarops/ADP/Products/

unset -v latest
for dir in $basedir/*
do
    for TLE in $dir/*TLE.txt
    do
         [[ $TLE -nt $latest ]] && latest=$TLE
    done
done
echo $latest
ln -s $latest Latest_TLE.txt



