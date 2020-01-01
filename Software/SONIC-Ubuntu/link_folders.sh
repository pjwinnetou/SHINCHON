INPUTDIR=$1
FOLDERSUFFIX=$2

if [ "$FOLDERSUFFIX" == "" ]
then
    FOLDERSUFFIX=$INPUTDIR
    echo $FOLDERSUFFIX
    echo
fi

for file in $INPUTDIR/*--*
do 
    centralityfolder=${file#*/}   
    echo $file $FOLDERSUFFIX-$centralityfolder
    ln -s $file $FOLDERSUFFIX-$centralityfolder
        
    cp master.params $FOLDERSUFFIX-$centralityfolder/
    cp master.trans $FOLDERSUFFIX-$centralityfolder/
done

