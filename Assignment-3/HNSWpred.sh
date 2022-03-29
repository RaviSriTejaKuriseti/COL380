#!/bin/sh

# Author : Ravi Sri Teja Kuriseti
# Script follows here:
mpirun --bind-to none ./wtf2 $1 $2 $3 $4
./wtf3
rm -r output.dat helper_n.txt helper.txt