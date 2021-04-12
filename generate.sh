#!/bin/bash

networkSize=$1
g=$2
ensembleSize=$3
coreNum=$4

name=N${networkSize}G${g}E${ensembleSize}C${coreNum}

g++ -O3 -march=native -flto -std=c++17 -o bin/${name}.out main-generate.cpp

./bin/${name}.out ${networkSize} ${g} ${ensembleSize} ${coreNum}
rm bin/${name}.out