#!/bin/bash

ColorOption=0
Directory="RunData"

if [ $1 ]
then 
    ColorOption=$1
fi

if [ $2 ] 
then
    Directory=$2
fi

./bin/vistools -colorOption=$ColorOption -directory=$Directory
python render.py -dir=$Directory