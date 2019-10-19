#!/bin/sh

#ulimit -s unlimited
ulimit -s 65500
mkdir run
cd run
mkdir output
ln -s ../build/2dvorticity_psm.exe ./
./2dvorticity_psm.exe
