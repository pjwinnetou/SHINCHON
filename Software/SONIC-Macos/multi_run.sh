#!/bin/bash
ROOTDIR=$(pwd)
TEMPLATEDIR="template"
INPUTDIR=$1
INPUTFILE="inited_pp_"
TMPDIR="tmp"
RESULTDIR=$ROOTDIR/results/$INPUTDIR

### overall batch logfile, nothing to do with SONIC runs
logfile=submission.log

echo
echo "-=-=-=-=-=-=-=-=-=-=-=-=-"
echo
echo -e "\e[1mThis starts the batch submission of several energy density profiles for the SONIC simulation package on the Eridanus Cluster.\e[0m"
echo
echo "-=-=-=-=-=-=-=-=-=-=-=-=-"
echo


#echo "Current directory: $ROOTDIR"

### Removes unnecessary slashes at the end of the CMD line argument
### Useful in conjunction with bash autocompletion
if [[ "${INPUTDIR:$((${#INPUTDIR}-1)):1}" == "/" ]]
then
    echo "Removing last character from CMD line argument '$INPUTDIR' because it is a slash '/'."
    echo
    INPUTDIR=${INPUTDIR%?}
    sleep 2
fi

#echo "Folder for initial energy densities: $INPUTDIR"
echo | tee -a $logfile
echo "$(date)" | tee -a $logfile
echo -e "Initial energy densities: \t $ROOTDIR/$INPUTDIR" | tee -a $logfile
echo -e "Current directory: \t\t $ROOTDIR" | tee -a $logfile
echo -e "Results directory: \t\t $RESULTDIR" | tee -a $logfile
echo

mkdir $RESULTDIR
status_resultdir=$?

cp $INPUTDIR/master.trans $RESULTDIR/
cp $INPUTDIR/master.params $RESULTDIR/

### checks whether $RESULTDIR is already present
if [[ $status_resultdir -eq 0 ]]
then

    echo -e "\e[1mSubmitting jobs.\e[0m"
    echo

    ### Extracting event number out of the file name
    ### Note the star '*' at the end of the folder
    for file  in $INPUTDIR/$INPUTFILE*
    do

	### Removes everything like '*_pp_', beginning of file name
	NUM=${file#*_pp_}
	### Removes everything like '.dat', end of file name
	NUM=${NUM%.dat}

	### Main submission routine of the jobs
	cp $TEMPLATEDIR/multi_vh2.pbs .
	sleep 1

	echo "Submitting initial condition: $INPUTFILE$NUM" | tee -a $logfile

	qsub -N VH2-${INPUTDIR}.$NUM -v TEMPLATEDIR=$TEMPLATEDIR,TMPDIR=$TMPDIR,INPUTDIR=$INPUTDIR,INPUTFILE=$INPUTFILE,NUM=$NUM,ROOTDIR=$ROOTDIR,RESULTDIR=$RESULTDIR multi_vh2.pbs | tee -a $logfile
	rm  multi_vh2.pbs

    done

    echo
    echo -e "\e[33mAll initial energy densities in '$INPUTDIR' submitted.\e[0m"
    echo "All initial energy densities in '$INPUTDIR' submitted." >> $logfile
    echo | tee -a $logfile

else 
    echo
    echo -e "\e[31mAborting:\e[0m Results folder already exists."
    echo "Aborting:\e[0m Results folder already exists." >> $logfile
    echo 
    exit
fi
