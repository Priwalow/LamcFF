#!/bin/bash

exp=$1
run=$2
run1=$3                                                                                                                                                                                               
(
check_process_url "http://bweb3/montecarlo.php?ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&ty=evtgen-charm&dt=Any&bl=caseB&dv=zfserv"
check_process_url "http://bweb3/montecarlo.php?ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&ty=evtgen-charged&dt=Any&bl=caseB&dv=zfserv"
check_process_url "http://bweb3/montecarlo.php?ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&ty=evtgen-mixed&dt=Any&bl=caseB&dv=zfserv"
)|wc -l 
