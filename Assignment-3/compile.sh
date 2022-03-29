#!/bin/sh

# Author : Ravi Sri Teja Kuriseti
# Script follows here:
g++ converter.cpp -o wtf1
mpic++ -std=c++17 main.cpp -o wtf2 -fopenmp
g++ fcon.cpp -o wtf3