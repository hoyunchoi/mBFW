#! /bin/bash

srcDir=/pds/pds11/hoyun/mBFW/src
libDir=/pds/pds11/hoyun/mBFW/lib
binDir=/pds/pds11/hoyun/mBFW/bin
library=/pds/pds11/hoyun/library

N=$1
G=$2
E=$3
C=$4

name=N${N}G${G}E${E}C${C}

function debugBuild {
	g++ -std=c++17 -Wall -g -fsanitize=address\
        -I ${libDir} -I ${library}\
        -o ${binDir}/${name}\
        ${srcDir}/main-generate.cpp
}

function build {
	g++ -std=c++17 -O3 -flto -march=native\
        -I ${libDir} -I ${library}\
        -o ${binDir}/${name}\
        ${srcDir}/main-generate.cpp
}

#* Compile
# debugBuild
build

#* Run
${binDir}/${name} $N $G $E $C
rm ${binDir}/${name}