### This runs VH2, after its completion B3D, after its completion Analyze                                                             
logdir=logdir
TEMPLATEDIR=template
RUNDIR=$1

if [[ $RUNDIR == "" ]]
then
    echo "No RUNDIR defined. Aborting."
    exit
fi

rm -rf $RUNDIR
cp -r $TEMPLATEDIR $RUNDIR

cd $RUNDIR

vh2id=$(qsub single_vh2.pbs)
echo $vh2id
#b3did=$(qsub -W depend=afterok:$vh2id mypbs_b3d)
#echo $b3did
#analyzeid=$(qsub -W depend=afterok:$b3did mypbs_analyze)
#echo $analyzeid
