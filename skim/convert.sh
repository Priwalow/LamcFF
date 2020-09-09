#!/bin/sh


if test $1
then
    maxproc=$1
else
    echo maximal amount of procecceses is not set 
    echo "                   " setting it to 50 
    maxproc=50
fi


for indexfile in `ls index/$2*`
do
    nproc=`ps -u privalov|wc -l`
    while  (( $nproc > $maxproc )) 
    do 
        sleep 1 
        nproc=`ps -u privalov|wc -l`
    done
	    file=`basename $indexfile .index`
	    echo running: $file   processes: $nproc / $maxproc
	    nice bfcp index/$file.index  mdst/$file.mdst >&mdst_log/$file.log & 
done

