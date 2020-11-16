#!/bin/bash

exp=$1
run=$2
run1=$3                                                                                                                                                                                               

check_process_url "http://bweb3/montecarlo.php?ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&ty=Any&dt=Any&bl=caseB&dv=zfserv" > mc_files.txt

while read i
do 
    echo `expr match "$i" '\</*.mdst\>'`
done < mc_files.txt
