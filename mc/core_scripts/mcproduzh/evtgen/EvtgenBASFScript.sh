#!/bin/bash

#Set env variables
export BELLE_LEVEL=b20090127_0910
export BELLE_DEBUG=opt
export CERN_LEVEL=2006
source /sw/belle/local/etc/bashrc_general

#Read input arguments
type=$1
decay_file=$2
evts=$3
Output_Name=$4

dat_file=${Output_Name}.mdst
log_file=${Output_Name}.log

#Run BASF
basf <<EOF >& ${log_file}

module register evtgen
path create main
path add_module main evtgen
module put_parameter evtgen USER_DECAY\\${decay_file}

$(cat $type)

initialize

table savebr gen_info
table savebr gen_beam
table save belle_event
table save mctype_all
table save hepevt_all

output open ${dat_file}
generate_event ${evts}
terminate
EOF
