#!/bin/bash

([ -z $1 ] || [ -z $2 ]) && echo "Usage: $0 filename num_parts" && exit

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

KAHYPAR=${SCRIPT_PATH}/cmake-build-release/kahypar/application/KaHyPar
CONFIG=${SCRIPT_PATH}/config/kaffpaD.ini
GRAPH=$1
K=$2
EPS=0.03

echo COMMAND: $KAHYPAR -h $GRAPH -k $K -e $EPS -o km1 -m acyclic -p $CONFIG --seed 0 --binary-kaffpaD=true
$KAHYPAR -h $GRAPH -k $K -e $EPS -o km1 -m acyclic -p $CONFIG --seed 0 --binary-kaffpaD=true
