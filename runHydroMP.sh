#!/bin/bash

keyInputs=$1

IFS="|" read -r var1 var2 var3 <<< "$keyInputs"

#export OMP_NUM_THREADS=$var1
#set OMP_NUM_THREADS $var1

folder_config="$var2"
folder_output="$var3"

#echo "Number Threads : ${OMP_NUM_THREADS}"

echo "Use config file in ${folder_config}/ and write results in ${folder_output}/ ..."

./beshydro --config ${folder_config} -o ${folder_output} -h