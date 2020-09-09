#!/bin/sh


(
for i in 07 09 11 13 15 17 19 21 23 25 27
do 
./run_exp.sh $i s
done
)|/bin/bash

