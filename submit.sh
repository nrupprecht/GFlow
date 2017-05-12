#!/bin/bash

LoadFile=HD-D20-Mu0.5-4x128-11.config

./driver -loadBuoyancy=$LoadFile -bR=0.2 -animate -writeDirectory=LongTube128 -CV -velocity=12 -time=10 -forces > longtube128.txt
