#!/bin/bash

exp=$1
run=$2
run1=$3                                                                                                                                                                                               


for i in $(check_process_url "http://bweb3/montecarlo.php?ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&ty=Any&dt=Any&bl=caseB&dv=zfserv")
do 
    echo "$i"
done 
