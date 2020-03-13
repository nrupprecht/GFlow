#!/bin/bash

ColorOption=0
Directory="RunData"
Data="Pos-1"

if [ $1 ]
then 
    ColorOption=$1
fi

if [ $2 ] 
then
    Directory=$2
fi

if [ $3 ]
then
    Data=$3
fi

./bin/vistools -colorOption=$ColorOption -directory=$Directory -data=$Data

python render.py -dir=$Directory -data=$Data
