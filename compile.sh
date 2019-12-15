#!/bin/bash
git submodule update --remote --init
mkdir cmake-build-release
cd cmake-build-release
cmake .. -DCMAKE_BUILD_TYPE=Release
make KaHyPar
cd ..