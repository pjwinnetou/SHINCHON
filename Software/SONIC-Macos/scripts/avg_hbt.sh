#!/bin/bash
adir=$1

outname=$adir/averages/av_pion_hbt.dat
allnum="$adir/data_*"

echo "Averaging over HBT radii: star_hbt_pion.dat"
echo "Current arguments start: ./avg_hbt.sh $adir"
echo "Output: $outname"
echo "Example: ./avg_vn_fluc.sh auau-40--50 pion"

declare -A ptarr;
declare -A v0arr;
declare -A v1arr;
declare -A v2arr;

filenum=0

for f in $allnum
do
    num=${f#*data_}
    #echo "$num"
    if [ "$num" = "ave" ]; then
       echo "skipping average first event "
    else
	file="$adir/data_$num/results/star_hbt_pion.dat"
	echo "Processing $file"
	a=-1  #line-counter
	while read line
	do
#	echo "Reading line $a : $line"
	    b=0  #array type counter
	    if [ "$a" -gt "-1" ]; then  #if it's not line -1, start filling in
		if [ "$filenum" = "0" ]; then
		    v0arr[$a]=0
		    v1arr[$a]=0
		    v2arr[$a]=0		
		fi
		for mydummy in $line
		do
		    ent=$(printf '%.10f' $mydummy)
		    if [ "$b" = "0" ]; then
			ptarr[$a]=$ent
		    fi
		    if [ "$b" = "1" ]; then           
			v0arr[$a]=$(echo "${v0arr[$a]}+$ent" | bc -l)
		    fi
		    if [ "$b" = "2" ]; then           
			v1arr[$a]=$(echo "${v1arr[$a]}+$ent" | bc -l)
		    fi
		    if [ "$b" = "3" ]; then
			v2arr[$a]=$(echo "${v2arr[$a]}+$ent" | bc -l)
		    fi
		    let b=b+1
		done
	    fi
	    let a=a+1
	    done<$file
	let filenum=filenum+1
    fi
done


declare -A v0bar;
declare -A v1bar;
declare -A v2bar;


for ((c=0;c<$a;c++))
do
    v0bar[$c]=$(echo "${v0arr[$c]}/$filenum" | bc -l)    
    v1bar[$c]=$(echo "${v1arr[$c]}/$filenum" | bc -l) 
    v2bar[$c]=$(echo "${v2arr[$c]}/$filenum" | bc -l) 
done

dummycounter=0

#second pass through data to get standard deviation
for f in $allnum
do
    num=${f#*data_}
    #echo "$num"
    if [ "$num" = "ave" ]; then
       echo "skipping average first event "
    else
	file="$adir/data_$num/results/star_hbt_pion.dat"
	echo "$file, second pass"
	a=-1  #line-counter
	while read line
	do
#	echo "Reading line $a : $line"
	    b=0  #array type counter
	    if [ "$a" -gt "-1" ]; then  #if it's not line -1, start filling in
		for mydummy in $line
		do
		    if [ "$dummycounter" = "0" ]; then
			v0arr[$a]=0
			v1arr[$a]=0
			v2arr[$a]=0		  
                    fi
		    ent=$(printf '%.10f' $mydummy)
		    if [ "$b" = "0" ]; then
			ptarr[$a]=$ent
		    fi
		    if [ "$b" = "1" ]; then           
			v0arr[$a]=$(echo "${v0arr[$a]}+($ent-${v0bar[$a]})*($ent-${v0bar[$a]})" | bc -l)
		    fi
		    if [ "$b" = "2" ]; then           
			v1arr[$a]=$(echo "${v1arr[$a]}+($ent-${v1bar[$a]})*($ent-${v1bar[$a]})" | bc -l)
		    fi
		    if [ "$b" = "3" ]; then
			v2arr[$a]=$(echo "${v2arr[$a]}+($ent-${v2bar[$a]})*($ent-${v2bar[$a]})" | bc -l)
		    fi
		    let b=b+1
		done
	    fi	
	    let a=a+1
	    done<$file
	let dummycounter=dummycounter+1
    fi
done

mkdir "$adir/averages"


printf "#PT[GeV]\t Rout\t s_Rout\t\t Rside\t s_Rside\t Rlong\t s_Rlong\n" > $outname

for ((c=0;c<$a;c++))
do
    dummyp=$(echo "${ptarr[$c]}/1000." | bc -l)
    dummy0=$(echo "sqrt(${v0arr[$c]})/$filenum" | bc -l)
    dummy1=$(echo "sqrt(${v1arr[$c]})/$filenum" | bc -l)
    dummy2=$(echo "sqrt(${v2arr[$c]})/$filenum" | bc -l)    
    printf "%2.2f\t\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t\n" $dummyp ${v0bar[$c]} $dummy0 ${v1bar[$c]} $dummy1 ${v2bar[$c]} $dummy2 >> $outname
done

echo "Averaged a total of $filenum files"

#echo "v0=${v0bar[0]}, $dummy"
