#!/bin/bash
TEMPLATEDIR=template
VH2DIR=uvh2-1
VH2BRANCH=b3d_ppbar
B3DDIR=b3d
B3DBRANCH=parallel
RESULTDIR=results
TMPDIR=tmp

### Creates template folder
mkdir -p $TEMPLATEDIR/output/default/cent0to5
mkdir -p $TEMPLATEDIR/analysis/default/cent0to5/details
mkdir -p $TEMPLATEDIR/progdata

cp $VH2DIR/generate $TEMPLATEDIR/
cp $VH2DIR/calcS $TEMPLATEDIR/
cp $VH2DIR/initE $TEMPLATEDIR/
cp $VH2DIR/vh2* $TEMPLATEDIR/

cp -r $VH2DIR/data $TEMPLATEDIR/
cp -r $VH2DIR/input $TEMPLATEDIR/
mkdir -p $TEMPLATEDIR/logdir
mkdir -p $TEMPLATEDIR/data/results

cp $B3DDIR/analyze $TEMPLATEDIR/
cp $B3DDIR/b3d $TEMPLATEDIR/
cp $B3DDIR/qualifiers.dat $TEMPLATEDIR/

cp -r $B3DDIR/parameters $TEMPLATEDIR/
cp -r $B3DDIR/install/progdata/madai/resinfo $TEMPLATEDIR/progdata/resinfo

cd $TEMPLATEDIR

ln -sf vh2* vh2
ln -sf analysis/default/cent0to5/details results_analysis

cd ../

### Creates other folders for batchsubmission
mkdir -p $TMPDIR
mkdir -p $RESULTDIR

echo
echo -e "\e[33mSONIC is now deployed.\e[0m"
echo
