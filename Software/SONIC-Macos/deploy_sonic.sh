#!/bin/bash
TEMPLATEDIR=template
VH2DIR=uvh2-1
VH2BRANCH=b3d_ppbar
B3DDIR=b3d
B3DBRANCH=parallel
RESULTDIR=results
TMPDIR=tmp
SCRIPTDIR=scripts
USESSHKEY=true

echo
echo "-=-=-=-=-=-=-=-=-=-=-=-=-"
echo
echo "Deploying SONIC by downloading and compiling the source files from Bitbucket.org and creating folder structure"
echo
echo "-=-=-=-=-=-=-=-=-=-=-=-=-"
echo

USERNAME=""

while true; do
	read -p "Are you using an SSH key with Bitbucket (y/n) (default n)? " RESP
  if [ -z "$RESP" ] || [ "$RESP" = "n" ] ; then
		USESSHKEY=false
		read -p "Enter your Bitbucket username: " USERNAME
		break
	elif [ "$RESP" = "y" ]
	then
		break
  else
		echo "Unrecognized answer. Please answer y or n."
	fi
done

### Downloads and compiles UVH2+1/VH2
if [ "$USESSHKEY" = true ]
then
    git clone -b $VH2BRANCH git@bitbucket.org:mluzum/uvh2-1.git $VH2DIR
else
    git clone -b $VH2BRANCH https://$USERNAME@bitbucket.org/mluzum/uvh2-1.git $VH2DIR
fi

cd $VH2DIR
make generate
make initE
make calcS
make vh2

cd ../

### Downloads and compiles B3D
if [ "$USESSHKEY" = true ]
then
    git clone -b $B3DBRANCH git@bitbucket.org:mhabich/b3d_sonic.git $B3DDIR
else
    git clone -b $B3DBRANCH https://$USERNAME@bitbucket.org/mhabich/b3d_sonic.git $B3DDIR
fi

cd $B3DDIR
make

cd ../

### Creates template folder
mkdir $TEMPLATEDIR
mkdir $TEMPLATEDIR/output
mkdir $TEMPLATEDIR/output/default
mkdir $TEMPLATEDIR/output/default/cent0to5

mkdir $TEMPLATEDIR/analysis
mkdir $TEMPLATEDIR/analysis/default
mkdir $TEMPLATEDIR/analysis/default/cent0to5
mkdir $TEMPLATEDIR/analysis/default/cent0to5/details

mkdir $TEMPLATEDIR/progdata

cp $VH2DIR/generate $TEMPLATEDIR/
cp $VH2DIR/calcS $TEMPLATEDIR/
cp $VH2DIR/initE $TEMPLATEDIR/
cp $VH2DIR/vh2* $TEMPLATEDIR/

cp -r $VH2DIR/data $TEMPLATEDIR/
cp -r $VH2DIR/input $TEMPLATEDIR/
mkdir $TEMPLATEDIR/logdir
mkdir $TEMPLATEDIR/data/results

cp $B3DDIR/analyze $TEMPLATEDIR/
cp $B3DDIR/b3d $TEMPLATEDIR/
cp $B3DDIR/qualifiers.dat $TEMPLATEDIR/

cp -r $B3DDIR/parameters $TEMPLATEDIR/
cp -r $B3DDIR/install/progdata/madai/resinfo $TEMPLATEDIR/progdata/resinfo

cd $TEMPLATEDIR

ln -s vh2* vh2
ln -s analysis/default/cent0to5/details results_analysis

cd ../

### Creates other folders for batchsubmission
mkdir $TMPDIR
mkdir $RESULTDIR

cp $SCRIPTDIR/single_run.sh .
cp $SCRIPTDIR/multi_run.sh .
cp $SCRIPTDIR/single_run.sh .
cp $SCRIPTDIR/link_folders.sh .
cp $SCRIPTDIR/single_*.pbs $TEMPLATEDIR/
cp $SCRIPTDIR/multi_*.pbs $TEMPLATEDIR/

cp $SCRIPTDIR/avg_*.sh $RESULTDIR
cp $SCRIPTDIR/int_qu.sh $RESULTDIR
cp $SCRIPTDIR/all_avg_all.sh $RESULTDIR
cp $SCRIPTDIR/all_int_qu.sh $RESULTDIR
cp $SCRIPTDIR/arrays.sh $RESULTDIR

echo
echo -e "\e[33mSONIC is now deployed.\e[0m"
echo
