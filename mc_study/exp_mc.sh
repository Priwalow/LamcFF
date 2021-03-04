#!/bin/sh


(


for i in 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do 
./run_exp_mc.sh $i l
done



for i in 07 09 11 13 15 17 19 21 23 25 27
do 
./run_exp_mc.sh $i s
done
)|/bin/bash
