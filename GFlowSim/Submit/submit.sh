#!/bin/bash

# Specify queue
# #$ -q long

# Specify job name
# #$ -N Stretch_Break_1000

# Number of iterations
iters=1000

# Set the label and tag
label="stretch_break"
tag=1
detection=0

# Force to apply
force=0.5

# Look for tag
if [ $1 ]
then
    tag=$1
fi

file="../Configurations/configuration_${label}_${tag}.config"
directory="../Configurations/data_${label}_${tag}"

# Create a random configuration
./../bin/gflow -config=../samples/pillar.cfg -phi=2 -time=120 -saveFile=$file -nowrite -quiet

# Run for some number of iterations
./../bin/gflow -time=120 -fracture -loadFile=$file -maxIters=$iters -dFy=$force -dFx=0 -writeDirectory=$directory -detection=$detection

# Clean up the configuration
rm $file
