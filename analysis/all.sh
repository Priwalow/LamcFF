#!/bin/sh


if test $1
then
    maxproc=$1
else
    echo maximal amount of procecceses is not set 
    echo "                   " setting it to 50 
    maxproc=50
fi


for mdstfile in `ls /gpfs/home/belle/uglov/LamcFF/mdst/$2*`
do
    nproc=`ps -u privalov|wc -l`
    while  (( $nproc > $maxproc )) 
    do 
	sleep 1 
	nproc=`ps -u privalov|wc -l`
    done

    
    file=`basename $mdstfile .mdst`
    echo running: $file   processes: $nproc / $maxproc
    if test -f /gpfs/home/belle/uglov/LamcFF/mdst/$file.mdst
	    then
	rm */$file.* >&/dev/null
	nice   ./run_mdst.sh $file &
    fi
done

