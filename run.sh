#!/bin/bash
type=$1
networkSize=$2
g=$3
ensembleSize=$4
machine=$5
coreNum=$6

spg run ${machine} ./compile/${type} ${networkSize} ${g} ${ensembleSize} ${machine} ${coreNum} >> ./log.txt