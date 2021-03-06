#!/bin/bash

srcDir=src
libDir=lib
binDir=bin
common=../library

networkSize=$1
g=$2
ensembleSize=$3
coreNum=$4

name=N${networkSize}G${g}E${ensembleSize}C${coreNum}

function debugBuild {
	g++ -std=c++17 -Wall -g -fsanitize=address\
        -I ${common} -I ${libDir}\
        -o ${binDir}/${name}\
	    ${srcDir}/main-generate.cpp
}

function build {
	g++ -std=c++17 -O3 -flto -march=native\
        -I ${common} -I ${libDir}\
        -o ${binDir}/${name} \
		${srcDir}/main-generate.cpp
}

#* Compile the source files
# build
debugBuild

#* Run
./${binDir}/${name} ${networkSize} ${g} ${ensembleSize} ${coreNum} >> log/generate.log
rm ${binDir}/${name}