FOLDERSUFFIX=$1

if [ "$FOLDERSUFFIX" == "" ]
then
    echo "No folder suffix defined. Exiting."
    exit
fi

for file in $FOLDERSUFFIX*
do 
    unlink $file
done

