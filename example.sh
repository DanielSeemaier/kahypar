#!/usr/bin/bash
git submodule update --remote --init
mkdir cmake-build-release
cd cmake-build-release
cmake .. -DCMAKE_BUILD_TYPE=Release
make kahypar
cd ..
gcc a.c -Iinclude/ -Lcmake-build-release/lib -l kahypar
LD_LIBRARY_PATH=cmake-build-release/lib
./a.out
