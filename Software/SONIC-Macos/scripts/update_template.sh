TEMPLATEDIR=template
VH2DIR=uvh2-1
VH2BRANCH=b3d_ppbar
B3DDIR=b3d
B3DBRANCH=parallel
SCRIPTDIR=scripts

cd $VH2DIR
make generate
make initE
make calcS
make vh2

cd ../

cd $B3DDIR
make

cd ../

cp $VH2DIR/generate $TEMPLATEDIR/
cp $VH2DIR/calcS $TEMPLATEDIR/
cp $VH2DIR/initE $TEMPLATEDIR/
cp $VH2DIR/vh2* $TEMPLATEDIR/

cp $B3DDIR/analyze $TEMPLATEDIR/
cp $B3DDIR/b3d $TEMPLATEDIR/
cp $B3DDIR/qualifiers.dat $TEMPLATEDIR/

echo
echo -e "\e[33mSONIC's executables have been updated.\e[0m"
echo