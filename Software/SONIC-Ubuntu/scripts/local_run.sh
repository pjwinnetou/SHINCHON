### This runs VH2, after its completion B3D, after its completion Analyze                                                             
logdir=logdir
TEMPLATEDIR=template
RUNDIR=$1

if [[ $RUNDIR == "" ]]
then
    echo "No RUNDIR defined. Aborting."
    exit
fi

cp -rn $TEMPLATEDIR $RUNDIR

cd $RUNDIR

echo "Running complete SONIC simulation."

./generate > $logdir/generate.log
./initE > $logdir/initE.log
./vh2 > $logdir/uvh2+1.log
./b3d default > $logdir/b3d.log
./analyze default > $logdir/b3d-analyze.log

