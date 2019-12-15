#!/bin/bash
git submodule update --remote --init
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make KaHyPar
cd ..