#!/bin/bash
runname=$1
echo "Trying to do cascade for run named $runname"
./b3d $runname
./averageall $runname
rm -fR ~/Dropbox/cascadeETA/$runname
mkdir ~/Dropbox/cascadeETA/$runname
cp output/$runname/cent0to5/meta.dat ~/Dropbox/cascadeETA/$runname/
cp output/$runname/cent0to5/av_Tabhist.dat ~/Dropbox/cascadeETA/$runname/
