#!/bin/bash

ColorOption=0

if [ $1 ]
then 
    ColorOption=$1
fi

./bin/vistools -colorOption=$ColorOption
python render.py