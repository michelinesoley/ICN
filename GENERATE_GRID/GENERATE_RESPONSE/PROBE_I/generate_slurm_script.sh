###################################
### Generate slurm input files	###
###################################

rm run_script_*.sh

num_files=4410
num_proc=10

files_proc=$(($num_files/$num_proc))

echo $files_proc

for proc in `seq 1 $num_proc`;
do
  scriptname=run_script_$proc.sh
  echo $scriptname

  offset=$((($proc-1)*$files_proc))

  cp slurm_template.inp $scriptname

  for line in `seq 1 $files_proc`;
  do
    
    label=$(($line+$offset-1))

    if [ $label -lt 10 ]
    then
      filename=qchem_000$label
    elif [ $label -lt 100 ]
    then
      filename=qchem_00$label
    elif [ $label -lt 1000 ]
    then
      filename=qchem_0$label
    else
      filename=qchem_$label
    fi

    echo $label
    echo $filename

    echo "qchem -nt \$SLURM_NTASKS $filename.in $filename.out" >> $scriptname

  done 

  echo " " >> $scriptname
  echo "echo 'Finishing program on'" >> $scriptname
  echo "date" >> $scriptname
  echo " " >> $scriptname

done

### change header

for file in run_script*
do

  echo $file

  prefix=$(basename $file .sh)

#  echo $prefix
  sed -i.bk "s/startFile/($prefix)_I/g" "$prefix"".sh"

done

rm *.bk

