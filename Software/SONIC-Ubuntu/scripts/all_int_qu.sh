#!/bin/bash
FOLDERPREFIX=$1

INTDIR=${FOLDERPREFIX}int

### Providing the array information
source arrays.sh

echo "Integrating over flow coefficients for ALL present centralities of a give collision system to calculate integrated quantities (dN/dy, <p_T>, etc.)"
echo "Current argument starts: ./all_int_qu.sh $FOLDERPREFIX"
echo "Output: $INTDIR"
echo "Example: ./all_int_qu.sh auau-"
echo 
echo

if [ "$FOLDERPREFIX" == "" ]
then
    echo "No folder suffix defined. Exiting."
    exit
fi

mkdir -p $INTDIR

### Writes headers for integrated quantities
for ((j=0;j<$N_species;j++))
do
    for ((i=0;i<$N_quantity;i++))
    do
	#outputfile="${quantity_file[$i]}_${species[$j]}.dat"
	#echo "$outputfile"
	outputfile="$INTDIR/${quantity_file[$i]}_${species[$j]}.dat"
	#echo -e  "#ctr_start \t ctr_mid \t ctr_end \t ${species[$j]}:${quantity_name[$i]}" > "${quantity_file[$i]}_${species[$j]}.dat"
	echo -e  "#ctr_start \t ctr_mid \t ctr_end \t ${species[$j]}:${quantity_name[$i]} \t stat_${quantity_name[$i]}"  > $outputfile
    done    
done


### Goes through all centrality classes of a particular system, labeled with the same FOLDERPREFIX
for file in $FOLDERPREFIX[0-9][0-9]--[0-9][0-9]{0,}
do

    centrality=${file#$FOLDERPREFIX}
    ### Loops over all species for a given centrality class
    for ((i=0;i<5;i++))
    do
	### $i corresponds to the species defined in 'arrays.sh'
	./int_qu.sh $FOLDERPREFIX $centrality $i
    done
done
