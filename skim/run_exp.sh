#!/bin/bash

if test $1
then
    exp=$1
else
    echo exp number is missing
    exit 1
fi



if test $2
then
    qtype=$2
else
    echo Queue type is not set, use \'s\' by default
    qtype=s
fi


for ((i=0; i<30; i++)) do
for ((j=0; j<10; j++)) do

if [ `./check_statistics.sh  $exp $i $j` -ne "0" ] 
then
echo bsub -q $qtype ./data.sh $exp ${i} ${j}
#bsub -q $qtype ./data.sh $exp ${i} ${j}
fi
done
done
