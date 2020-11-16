#!/bin/bash -f

BL=$1
exp=$2
basf=$3
expno6=$4
run=$5
input=$6
outfile=$7

set -x

export BELLE_LEVEL=$BL
export BELLE_DEBUG=opt
if [ ${expno6#0} -ge 31 ]
then
	export CERN_LEVEL=2006
fi
export USE_GRAND_REPROCESS_DATA=1
set +x
source /sw/belle/local/etc/bashrc_general
set -x
export BASF_NPROCESS=4
export ADDBG_EXP=$exp
export ADDBG_EXP6=$expno6
export ADDBG_RUN=$run
source .addbgrc

basf <<EOF >& ${outfile}.log
$(cat $basf)

histogram define ${outfile}.hbk
output open ${outfile}.mdst
process_event ${input}
terminate
EOF
