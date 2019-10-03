#!/bin/bash

exp=$1
run=$2
run1=$3                                                                                                                                                                                               
(
check_process_url "http://bweb3/mdst.php?bl=caseB&skm=HadronB&ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&dv=zfserv&dt=Any"
check_process_url "http://bweb3/mdst.php?bl=caseB&skm=HadronBJ&ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&dv=zfserv&dt=Any"
)|wc -l 
