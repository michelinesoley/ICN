###################################################################
### Generate input files for QCHEM for a sets of configurations	###
###################################################################

for file in ../../*.coor
do
  echo $file
  prefix=$(basename $file .coor)
  output=$prefix'.in'

  echo '$molecule' > $output
  cat $file >> $output
  echo '$end' >> $output
  cat input-trimed.in >> $output
done




