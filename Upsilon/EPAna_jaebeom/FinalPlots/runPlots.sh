
nRun='10'
kInitPos='1' # 1 : Upsilon initial position folllowing Initial Gaussian distribution
fitorder='3'

for ((i=0;i<2;i++));
do
  echo isLine $i
  root -l -q -b 'drawRAA_nPart_compare_data.C('$nRun','$kInitPos','$i')'
  root -l -q -b 'drawRAA_pt_compare_data.C('$nRun','$kInitPos','$i')'
  root -l -q -b 'drawV2_pt_compare_data.C('$fitorder','$nRun','$kInitPos','$i')'
  root -l -q -b 'drawV2_pt_compare_data.C('$fitorder','$nRun','$kInitPos','$i')'
  for ((id=0;id<3;id++));
  do
    echo data $id
    root -l -q -b 'drawRAA_nPart.C('$nRun','$kInitPos','$id','$i')'
    root -l -q -b 'drawV2_pt.C('$fitorder','$nRun','$kInitPos','$id','$i')'
    root -l -q -b 'drawV2_pt_compare.C('$fitorder','$nRun','$kInitPos','$id','$i')'
    root -l -q -b 'drawRAA_pt.C('$nRun','$kInitPos','$id','$i')'
    root -l -q -b 'drawRAA_pt.C('$nRun','$kInitPos','$id','$i')'
  done
done

