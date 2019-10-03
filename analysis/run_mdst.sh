#!/bin/sh

if test $1
then
    file=$1
else
    echo exp number is missing
    exit 1
fi

source /sw/belle/local/etc/bashrc_general



export USE_GRAND_REPROCESS_DATA=1
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so



#export BELLE_MESSAGE_LEVEL=DDEBUG
export BELLE_MESSAGE_LEVEL=INFO
#### case B

unset BELLE_USE_TMP

export QQ_USER_TABLE="./user.dec"


filelist=`find /gpfs/home/belle/uglov/LamCFF/mdst -name $file\*mdst -size +0c|sort`


outfile=hbk/$file.h
logfile=log/$file.log

(
cat <<EOF


path add_module main fix_mdst User_reco
path add_condition main <:0:KILL

  

initialize


histogram define $outfile




EOF

for file in ${filelist}
do
    echo "process_event ${file} " 
done


echo "terminate"
) |basf> $logfile 2>&1
#gzip -9 $outfile 