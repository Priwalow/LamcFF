if test $1
then
    exp=$1
else
    echo exp number is missing
    exit 1
fi


for ((i=0; i<30; i++)) do
for ((j=0; j<10; j++)) do

if [ `./check_MC.sh  $exp $i $j` -ne "0" ] 
then
echo $exp ${i} ${j}
fi
done
done
 
