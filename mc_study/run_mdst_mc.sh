#!/bin/sh

if test $1
then
    exp=$1
else
    echo exp number is missing
    exit 1
fi

if test $2
then
    run=$2
    else
    echo run number is missing
    exit 1
fi

if test $3
then
    run1=$3
    else
    echo run1 number is missing
    exit 1
fi


#####################################
source /sw/belle/local/etc/bashrc_general



export USE_GRAND_REPROCESS_DATA=1
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so



#export BELLE_MESSAGE_LEVEL=DDEBUG
export BELLE_MESSAGE_LEVEL=INFO
#### case B

unset BELLE_USE_TMP

export QQ_USER_TABLE="./user.dec"


(
cat <<EOF


path add_module main fix_mdst User_reco
path add_condition main <:0:KILL

  

initialize


histogram define hbk/${exp}.${run}.${run1}.h

process_url "http://bweb3/montecarlo.php?ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&ty=evtgen-charm&dt=Any&bl=caseB&dv=zfserv"



EOF

echo terminate

) |basf >  ./log/${exp}.${run}.${run1} 2>&1
#gzip -9 hbk/${exp}.${run}.${run1}.h
