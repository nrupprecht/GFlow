#!/bin/bash

num=200
if [ $1 ]
then
    num=$1
fi

# Get the starting 
start=1
if [ $2 ]
then 
    start=$2
fi

for i in `seq $start $num`;
do
    qsub submit.sh $i
done
