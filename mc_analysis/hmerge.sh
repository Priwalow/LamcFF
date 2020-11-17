#!/bin/sh

echo "#!/bin/sh" > zzzz_hmerge.sh
#echo "cd hmerge" >> zzzz_hmerge.sh

for i in  07 09  11 13 15 17 19 21 23 25 27  31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do

(
echo hmerge hmerge/$i.h _
ls hbk/$i.*.h |awk '{print $1,"_"}'
echo 

) >hmerge/$i.kumac
echo paw -b hmerge/$i.kumac >> zzzz_hmerge.sh "&"

done


sh zzzz_hmerge.sh
sleep 3 
#rm zzzz_hmerge.sh
#rm hmerge/*kumac