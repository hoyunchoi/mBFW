#!/bin/bash

networkSize=$1
g=$2

name=N${networkSize}G${g}

g++ -O3 -march=native -flto -std=c++17 -o bin/${name} main-data_class.cpp

./bin/${name} ${networkSize} ${g} > data.log
rm bin/${name}
