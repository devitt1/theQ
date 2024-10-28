#!/bin/bash
echo "theQ with Quplexity integration."
as ./Quplexity/ARM/math.s -o ./Quplexity/ARM/math.o
gcc main.c Simulator/sim.c Simulator/norm.c -lm Quplexity/ARM/math.o -o sim
./sim