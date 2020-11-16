#!/bin/tcsh -f

if( $2 == '' ) then
  echo "Usage : $0 [gen] [expno] {runno}"
  exit
endif

set genfile = $1
if( -e $genfile ) then
else
  echo "$genfile not found"
  exit
endif

set exp = $2                                     # 7, 07, 19a etc.
set expno = `echo $exp | sed s/'[a-z]'//`        # 7, 07 19 etc.
set expdiv = `echo $exp | sed s/'[0-9]'//g`      #        a etc.
set expno2 = `printf "%02d" $expno` # 07 07 19 etc.
set expno6 = `printf "%06d" $expno` # 000007 000007 000019 etc
set basfexp = $expno2$expdiv
set gsimexp = $expno

if( $3 == '' ) then
  set run = `echo 110 | awk '{srand(); print int(rand()*$1)}'`
else
  set run = $3
endif

set bl = "b20060529_2127"
if( $expno <= 27 ) then
  set bl = "b20030807_1600"
else if( $expno <= 71 ) then
  set bl = "b20090127_0910"
else
  echo "Exp $expno is not implemented in rec.csh"
  exit
endif
#unsetenv MY_TOP_DIR
setenv BELLE_LEVEL $bl
setenv BELLE_DEBUG opt
if( $expno >= 31 ) then
 setenv CERN_LEVEL 2006
endif
setenv USE_GRAND_REPROCESS_DATA 1
setenv BASF_NPROCESS 4
source /sw/belle/local/etc/cshrc_general
#setenv BELLE_INDEX_POSTGRES_MC_SERVER bsql3
#setenv BELLE_MESSAGE_LEVEL INFO

set basf = "basf/mcnobg-e$basfexp.basf"
set gsim = "gsim/gsim.$gsimexp.dat"

if( !(-e $basf) ) then
  echo "$basf not found"
  exit
endif
if( !(-e $gsim) ) then
  echo "$gsim not found"
  exit
endif

set logdir = "./log"
mkdir -p $logdir
set logfile = $logdir/$genfile:t:r.log

set hbkdir = "./hbook"
mkdir -p $hbkdir
set hbkfile = $hbkdir/$genfile:t:r.hbook

set outdir = "./mdst"
mkdir -p $outdir
set outfile = $outdir/$genfile:t:r.mdst

#setenv ADDBG_EXP  $exp
#setenv ADDBG_EXP6 $expno6
#setenv ADDBG_RUN  $run
#source .addbgrc

#echo run = $run
#echo hbk = $hbkfile
#echo out = $outfile
#echo gen = $genfile
#echo $BHOME
#echo $BELLE_LEVEL
#ls $basf
#ls $gsim
#which basf
#ls  $ADDBG_DAT
#echo $BELLE_RUN_DIR
#exit

basf <<EOF >& $logfile
`cat $basf`
histogram define $hbkfile
output open $outfile
process_event $genfile
terminate
EOF
