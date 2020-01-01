#!/bin/bash
adir=$1
species=$2
cumulants=$3

outname=$(printf '%s/averages/av_%s_vn-%s.dat' $adir $species $cumulants)
allnum="$adir/data_*"

echo "Averaging over flow coefficients stemming from cumulants: ..._vn_qX.dat"
echo "Current arguments start: ./avg_vn_q.sh $adir $species $cumulants"
echo "Output: $outname"
echo "Example: ./avg_vn_q.sh auau-40--50 pion 4th"

declare -A ptarr;
declare -A v0arr;
declare -A v1arr;
declare -A v2arr;
declare -A v3arr;
declare -A v4arr;
declare -A v5arr;
declare -A v6arr;

filenum=0

for f in $allnum
do
    num=${f#*data_}
    #num=${tmp%.*}
    echo "$num"
    if [ "$num" = "ave" ]; then
       echo "skipping average first event "
    else
	file=$(printf '%s/data_%s/results/%s_vn-%s.dat' $adir $num $species $cumulants)
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
		    v3arr[$a]=0
		    v4arr[$a]=0
		    v5arr[$a]=0
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
		    if [ "$b" = "4" ]; then
			v3arr[$a]=$(echo "${v3arr[$a]}+$ent" | bc -l)
		    fi
		    if [ "$b" = "5" ]; then
			v4arr[$a]=$(echo "${v4arr[$a]}+$ent" | bc -l)
		    fi
		    if [ "$b" = "6" ]; then
			v5arr[$a]=$(echo "${v5arr[$a]}+$ent" | bc -l)
		    fi
#		if [ "$b" = "6" ]; then
#		    v6arr[$a]=$(echo "${v6arr[$a]}+$ent" | bc -l)
#		fi
		    let b=b+1
		done
	    fi
#	if [ "$a" = "6" ]; then
#          echo "${ptarr[$a]}, v0=${v0arr[$a]} v1=${v1arr[$a]} v2=${v2arr[$a]}"
#	fi
	    let a=a+1
	    done<$file
	let filenum=filenum+1
    fi
done

#echo "total $filenum, v1=${v1arr[0]}"



declare -A v0bar;
declare -A v1bar;
declare -A v2bar;
declare -A v3bar;
declare -A v4bar;
declare -A v5bar;

for ((c=0;c<$a;c++))
do
    v0bar[$c]=$(echo "${v0arr[$c]}/$filenum" | bc -l)    
    v1bar[$c]=$(echo "${v1arr[$c]}/$filenum" | bc -l) 
    v2bar[$c]=$(echo "${v2arr[$c]}/$filenum" | bc -l)
    v3bar[$c]=$(echo "${v3arr[$c]}/$filenum" | bc -l) 
    v4bar[$c]=$(echo "${v4arr[$c]}/$filenum" | bc -l)
    v5bar[$c]=$(echo "${v5arr[$c]}/$filenum" | bc -l)
#    echo "${ptarr[$c]}  ${v0bar[$c]}  ${v1bar[$c]}   ${v2bar[$c]}  ${v3bar[$c]}  ${v4bar[$c]}   ${v5bar[$c]} " >> $outname
done

dummycounter=0

#second pass through data to get standard deviation
for f in $allnum
do
    num=${f#*data_}
    #num=${tmp%.*}
    #echo "$num"
    if [ "$num" = "ave" ]; then
       echo "skipping average first event "
    else
	file=$(printf '%s/data_%s/results/%s_vn-%s.dat' $adir $num $species $cumulants)
    #file="$adir/data_$num/results/pion_vn_fluc.dat"
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
			v3arr[$a]=0
			v4arr[$a]=0
			v5arr[$a]=0
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
		    if [ "$b" = "4" ]; then
			v3arr[$a]=$(echo "${v3arr[$a]}+($ent-${v3bar[$a]})*($ent-${v3bar[$a]})" | bc -l)
		    fi
		    if [ "$b" = "5" ]; then
			v4arr[$a]=$(echo "${v4arr[$a]}+($ent-${v4bar[$a]})*($ent-${v4bar[$a]})" | bc -l)
		    fi
		    if [ "$b" = "6" ]; then
			v5arr[$a]=$(echo "${v5arr[$a]}+($ent-${v5bar[$a]})*($ent-${v5bar[$a]})" | bc -l)
		    fi
#		if [ "$b" = "6" ]; then
#		    v6arr[$a]=$(echo "${v6arr[$a]}+$ent" | bc -l)
#		fi
		    let b=b+1
		done
	    fi	
	    let a=a+1
	    done<$file
	let dummycounter=dummycounter+1
    fi
done

mkdir "$adir/averages"


printf "#PT[GeV]\t V0\t s_V0\t\t V1\t s_V1\t V2\t s_V2\t V3\t s_V3\t V4\t s_V4\t V5\t s_V5\n" > $outname

#echo "#  PT [GeV]           V0                       V1                       V2                       V3                        V4                     V5" > $outname


#printf "%f %10s %s\n" 5.5 '' [UP] >> $outname

for ((c=0;c<$a;c++))
do
    dummy0=$(echo "sqrt(${v0arr[$c]})/$filenum" | bc -l)
    dummy1=$(echo "sqrt(${v1arr[$c]})/$filenum" | bc -l)
    dummy2=$(echo "sqrt(${v2arr[$c]})/$filenum" | bc -l)
    dummy3=$(echo "sqrt(${v3arr[$c]})/$filenum" | bc -l)
    dummy4=$(echo "sqrt(${v4arr[$c]})/$filenum" | bc -l)
    dummy5=$(echo "sqrt(${v5arr[$c]})/$filenum" | bc -l)
    printf "%2.2f\t %10.4f\t %7.8f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\n" ${ptarr[$c]} ${v0bar[$c]} $dummy0 ${v1bar[$c]} $dummy1 ${v2bar[$c]} $dummy2 ${v3bar[$c]} $dummy3 ${v4bar[$c]} $dummy4 ${v5bar[$c]} $dummy5 >> $outname
done

echo "Averaged a total of $filenum files"

#echo "v0=${v0bar[0]}, $dummy"
