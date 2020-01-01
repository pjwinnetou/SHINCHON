execname=vh2-2.1
rm $execname.e*
rm $execname.o*
make clean
make initE
make vh2
./initE
echo "
--------------------------------------------------------------------------------------------------------
Completed:
make clean
make initE
make vh2
./initE"
echo "Running $execname in"
for (( time=5; time>=1;  time-- ))
do
echo "$time sec"
sleep 1
done
rm data/snapshot/*
qsub mypbs
#./$execname
