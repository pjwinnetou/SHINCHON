#!/bin/bash

FOLDERPREFIX=$1
centrality=$2
species_number=$3

### Providing the array information
source arrays.sh

species=${species[$species_number]}

if [ "$FOLDERPREFIX" == "" ]
then
    echo "No folder suffix defined. Exiting."
    exit
fi

### Actual folder on which script acts
adir=$FOLDERPREFIX$centrality

### Directory for the integrated quantities/results
INTDIR=${FOLDERPREFIX}int

### Interval of centrality class
START=${centrality%--*}
END=${centrality#*--}
MID=$(echo "scale=1; ($START+$END)/2" | bc )

echo "Integrating over flow coefficients to calculate integrated quantities (dN/dy, <p_T>, etc.)"
echo "Current argument starts: ./int_qu.sh $FOLDERPREFIX $centrality $species_number"
echo "Output: $INTDIR/int_..._$species.dat"
echo "Example: ./int_qu.sh auau- 40--50 2"
echo 
dndy=0;
tdndy=0;
cdndy=0;
sdndy=0;
mpt=0;
smpt=0;
mpt2=0;
dpt=0;

v2=0;
tv2=0;
sv2=0;
cv2=0;
scv2=0;
v3=0;
sv3=0;
v4=0;
sv4=0;

file=$(printf '%s/averages/av_%s_vn.dat' $adir $species)
echo "Processing $file"
a=-1  #line-counter
while read line
do
    b=0  #array type counter
    if [ "$a" -gt "-1" ]; then  #if it's not line -1, start filling in
	for mydummy in $line
	do
	    ent=$(printf '%.10f' $mydummy)	   
	    if [ "$b" = "0" ]; then
		if [ "$a" = "1" ]; then
		    dpt=$(echo "$ent-$pt" | bc -l)
		fi
		pt=$ent		
	    elif [ "$b" = "1" ]; then
		v0=$ent
		#echo "v0=$v0"
		dndy=$(echo "$dndy+$v0*$pt" | bc -l)
		if [ "$a" = "0" ]; then
		    dndy=$(echo "$dndy-$v0*$pt*0.5" | bc -l)
		fi		
		mpt=$(echo "$mpt+$v0*$pt*$pt" | bc -l)
		mpt2=$(echo "$mpt2+$v0*$pt*$pt*$pt" | bc -l)
		if [ "$a" -gt "2" ]; then
                cdndy=$(echo "$cdndy+$v0*$pt" | bc -l)
		fi
		if [ "$a" = "3" ]; then
                    cdndy=$(echo "$cdndy-$v0*$pt*0.5" | bc -l)
                fi

	    elif [ "$b" = "2" ]; then
		sv0=$ent
		#echo "v0=$v0"
		sdndy=$(echo "$sdndy+$sv0*$pt" | bc -l)
		smpt=$(echo "$smpt+$sv0*$pt*$pt" | bc -l)		   
		if [ "$a" = "0" ]; then
                    sdndy=$(echo "$sdndy-$sv0*$pt*0.5" | bc -l)
		    smpt=$(echo "$smpt-$sv0*$pt*$pt*0.5" | bc -l)
                fi
	    elif [ "$b" = "5" ]; then
		v2=$(echo "$v2+$v0*$pt*$ent" | bc -l)
		if [ "$a" = "0" ]; then
		    v2=$(echo "$v2-$v0*$pt*$ent*0.5" | bc -l)
		fi
		if [ "$a" -gt "2" ]; then
                cv2=$(echo "$cv2+$v0*$pt*$ent" | bc -l)
		fi 
		if [ "$a" = "3" ]; then
                cv2=$(echo "$cv2-0.5*$v0*$pt*$ent" | bc -l)
                fi
	    elif [ "$b" = "6" ]; then		
		sv2=$(echo "$sv2+$v0*$pt*$ent" | bc -l)
		if [ "$a" = "0" ]; then
                    sv2=$(echo "$sv2-$v0*$pt*$ent*0.5" | bc -l)
                fi
	    elif [ "$b" = "7" ]; then
		v3=$(echo "$v3+$v0*$pt*$ent" | bc -l)
		if [ "$a" = "0" ]; then
                    v3=$(echo "$v3-$v0*$pt*$ent*0.5" | bc -l)
                fi
	    elif [ "$b" = "8" ]; then
		sv3=$(echo "$sv3+$v0*$pt*$ent" | bc -l)
		if [ "$a" = "0" ]; then
                    sv3=$(echo "$sv3-$v0*$pt*$ent*0.5" | bc -l)
                fi
	    elif [ "$b" = "9" ]; then
		v4=$(echo "$v4+$v0*$pt*$ent" | bc -l)
	    elif [ "$b" = "10" ]; then
		sv4=$(echo "$sv4+$v0*$pt*$ent" | bc -l)
	    fi
	    let b=b+1	    
	done
    fi
    let a=a+1
    done<$file


quantity[0]=$(echo "$dndy*$dpt*12.5664" | bc -l)
quantity_error[0]=$(echo "$sdndy*$dpt*12.5664" | bc -l)

quantity[1]=$(echo "$mpt/$dndy" | bc -l)
quantity_error[1]=$(echo "$smpt/$dndy" | bc -l)


quantity[2]=$(echo "$v2/$dndy" | bc -l)
quantity_error[2]=$(echo "$sv2/$dndy" | bc -l)

quantity[3]=$(echo "$cv2/$cdndy" | bc -l)

quantity[4]=$(echo "$v3/$dndy" | bc -l)
quantity_error[4]=$(echo "$sv3/$dndy" | bc -l)



### This need to be redefined in forms of the 'quantity' arrays, currently not output
#finalmpt2=$(echo "$mpt2/$dndy" | bc -l)
#finaltv2=$(echo "($v2+$tv2)/($dndy+$tdndy)" | bc -l)
#finaltdndy=$(echo "($dndy+$tdndy)*$dpt*12.5664" | bc -l)
#echo "where v2=$v2 tv2=$tv2" 
#echo "Total v2(trapezoid)=$finaltv2"

### Output routines

for ((i=0;i<$N_quantity;i++))
do
    outputfile="$INTDIR/${quantity_file[$i]}_${species[$j]}.dat"
    echo -e  "$START \t $MID  \t $END \t ${quantity[$i]} \t ${quantity_error[$i]}" >> $outputfile
done
