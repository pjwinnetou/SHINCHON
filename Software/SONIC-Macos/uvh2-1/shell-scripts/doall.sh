#!/bin/bash
runname=$1
echo "Trying to do all for run named $runname"
rm data/snapshot/*
rm data/results*
./initE
./vh2-1.3
cp -r parameters/default parameters/$runname
cp -r output/default output/$runname
cp -r data output/$runname/cent0to5/
cp data/freezeout.dat output/$runname/cent0to5/
cp data/meta.dat output/$runname/cent0to5/
cp etaoverS.dat output/$runname/cent0to5/
cp beta2.dat output/$runname/cent0to5/
cp lambda1.dat output/$runname/cent0to5/
./b3d $runname
./averageall $runname
rm -fR ~/Dropbox/cascadeETA/$runname
mkdir ~/Dropbox/cascadeETA/$runname
cp output/$runname/cent0to5/meta.dat ~/Dropbox/cascadeETA/$runname/
cp output/$runname/cent0to5/av_Tabhist.dat ~/Dropbox/cascadeETA/$runname/
