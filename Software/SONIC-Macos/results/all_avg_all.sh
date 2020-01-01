#!/bin/bash
FOLDERPREFIX=$1

### Providing the shared variables and arrays
source arrays.sh

if [ "$FOLDERPREFIX" == "" ]
then
    echo "No folder suffix defined. Exiting."
    exit
fi

echo "Averaging over HBT radii and flow coefficients stemming from fluctuations for ALL present centralities of a given collision system."
echo "Current arguments start: ./all_avg_all.sh $FOLDERPREFIX"
echo "Output: ${FOLDERPREFIX}xx--yy/averages"
echo "Example: ./all_avg_all.sh auau-"
echo "(Remove the centrality class (xx--yy) from the string"

### Goes through all centrality classes of a particular system, labeled with the same FOLDERPREFIX
for file in $FOLDERPREFIX[0-9][0-9]--[0-9][0-9]{0,}
do
    echo "./avg_all.sh $file"
    ./avg_all.sh $file
done

