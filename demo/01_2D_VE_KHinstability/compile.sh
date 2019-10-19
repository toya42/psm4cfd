#!/bin/bash
mkdir build
cd build
#cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=debug ../
cmake -D CMAKE_Fortran_COMPILER=ifort -D CMAKE_BUILD_TYPE=debug ../
make

