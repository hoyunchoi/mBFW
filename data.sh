#!/bin/bash

srcDir=src
libDir=lib
binDir=bin
common=../library

networkSize=$1
g=$2

name=N${networkSize}G${g}

function debugBuild {
	g++ -std=c++17 -Wall -g -fsanitize=leak -I ${common} \
		${srcDir}/main-data.cpp \
		-o ${binDir}/${name}
}

function build {
	g++ -std=c++17 -O3 -flto -march=native -I ${common} \
		${srcDir}/main-data.cpp \
		-o ${binDir}/${name}
}

#* Compile the source files
build
# debugBuild

#* Run
cd ${binDir}
./${name} ${networkSize} ${g} > ../log/data.log
rm ${name}
