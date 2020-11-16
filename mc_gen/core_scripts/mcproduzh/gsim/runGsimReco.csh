#! /bin/tcsh -f

# Arguments
# 1: absolute path to evtgen mdst files

if( $1 == '' ) then
  echo "Usage : ./runGsimReco.csh [absolutePathToEvtgenGeneratorFiles/]"
  exit
endif

set evtgenmdstpath=$1

echo "Will simulate detector response for all files in directory:"
echo $evtgenmdstpath

foreach file (`\ls -1 $evtgenmdstpath`)
   set exp=`echo $file | awk -F '_' '{print $3}'`

   set expno = `echo $exp | sed s/'0'//g`     

   echo $evtgenmdstpath$file $expno
   #bsub -q b_a ./rec.csh $evtgenmdstpath$file $expno
end

