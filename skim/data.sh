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

module register fix_mdst User_reco user_index
path add_module main fix_mdst  User_reco user_index
path add_condition main <:0:KILL

initialize
output open       ../index/${exp}.${run}.${run1}.index
histogram define  ../hbk/${exp}.${run}.${run1}.hist

process_url "http://bweb3/mdst.php?bl=caseB&skm=HadronB&ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&dv=zfserv&dt=Any"

process_url "http://bweb3/mdst.php?bl=caseB&skm=HadronBJ&ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&dv=zfserv&dt=Any"

EOF

echo terminate

) |basf >  ../log/${exp}.${run}.${run1} 2>&1

