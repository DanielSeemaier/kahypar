#!/bin/bash

# get directory of this shellscript
pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}";
while([ -h "${SCRIPT_PATH}" ]); do
    cd "`dirname "${SCRIPT_PATH}"`"
    SCRIPT_PATH="$(readlink "`basename "${SCRIPT_PATH}"`")";
done
cd "`dirname "${SCRIPT_PATH}"`" > /dev/null
SCRIPT_PATH="`pwd`";
popd  > /dev/null

KAHYPAR=${SCRIPT_PATH}/build/kahypar/application/KaHyPar
CONFIG=${SCRIPT_PATH}/config/kaffpaD.ini
GRAPH=$1
K=$2
EPS=0.03

$KAHYPAR -h $GRAPH -k $K -e $EPS -o km1 -m acyclic -p $CONFIG --seed 0 --shm