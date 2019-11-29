#!/bin/sh

ulimit -s unlimited
mkdir run
cd run
mkdir output
ln -s ../build/2dvorticity_psm.exe ./
./2dvorticity_psm.exe
