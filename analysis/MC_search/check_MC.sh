#!/bin/bash

exp=$1
run=$2
run1=$3                                                                                                                                                                                               
(
check_process_url "http://bweb3/montecarlo.php?bl=caseB&ex=${exp}$&rs=${run}${run1}0&re=${run}${run1}9&ty=Any&dv=zfserv&dt=Any"
)|wc -l 
