#! /bin/bash

srcDir=/pds/pds11/hoyun/mBFW/src
libDir=/pds/pds11/hoyun/mBFW/lib
binDir=/pds/pds11/hoyun/mBFW/bin
library=/pds/pds11/hoyun/library

name=test

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
${binDir}/${name}
rm ${binDir}/${name}