#!/bin/bash
### Directory on which script acts
RESULTDIR=$1

### Providing the shared variables and arrays
source arrays.sh

echo "Averaging over HBT radii and flow coefficients stemming from fluctuations: star_hbt_pion.dat, ..._vn_fluc.dat"

### Removes eventual slash '/' 
if [[ "${RESULTDIR:$((${#RESULTDIR}-1)):1}" == "/" ]]
then
    echo
    echo "Removing last character from CMD line argument '$RESULTDIR' because it is a slash '/'."
    echo
    RESULTDIR=${RESULTDIR%?}
    sleep 1
fi

echo "Current arguments start: ./avg_all.sh $RESULTDIR"
echo "Output: $RESULTDIR/averages"
echo "Example: ./avg_all.sh auau-40--50"

### Actual averaging procedure
for ((i=0;i<5;i++))
do
    ./avg_vn_fluc.sh $RESULTDIR ${species[$i]} &
done
./avg_hbt.sh $RESULTDIR 

### Experimental part
#./avg_vn_q.sh $RESULDIR