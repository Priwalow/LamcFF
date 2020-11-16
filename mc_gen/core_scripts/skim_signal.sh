#!/bin/bash
# Arguments
# 1: Name of analysis module
# 2: Absolute path to .mdst file
# 3: Output file (minus extension)

MODULE=$1
inputfile=$2
outputfile=$3

export BELLE_LEVEL=b20090127_0910
export BELLE_DEBUG=opt
export CERN_LEVEL=2006
source /sw/belle/local/etc/bashrc_general
export BASF_NPROCESS=0

basf <<EOF >& ${outputfile}.log
module register fix_mdst
module register ${MODULE}
path create main
path create analysis
path add_module main fix_mdst
path add_module main ${MODULE}
path add_condition main >:0:analysis
path add_condition main <=:0:KILL

initialize
table save mdst_all
table save evtcls_all
table save evtvtx_all
table save gsim_rand
table save hepevt_all
table save mctype
table save level4_all

table save bgtbl_info
table save dattof_trgl0
table save reccdc_timing

histogram define ${outputfile}.hbk
output open ${outputfile}.mdst
process_event ${inputfile} 0 
output close
terminate
EOF
	rm /group/belle/users/tbloom/MC_data/${folder}/${mdst}_skim.hbk
done
