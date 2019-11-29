#!/bin/bash

# 01 debug+OFF
mkdir build
cd build
cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=debug -D large_array=OFF ../
make
cd ..

ulimit -s unlimited
mkdir run
cd run
mkdir output
ln -s ../build/2dvorticity_psm.exe ./
./2dvorticity_psm.exe

# 02 debug+ON
cd ../build
make clean
cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=debug -D large_array=ON ../
cd ../run
./2dvorticity_psm.exe

# 03 fast+OFF
cd ../build
make clean
cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=fast -D large_array=OFF ../
cd ../run
./2dvorticity_psm.exe

# 03 fast+ON
cd ../build
make clean
cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=fast -D large_array=ON ../
cd ../run
./2dvorticity_psm.exe
cd ..
