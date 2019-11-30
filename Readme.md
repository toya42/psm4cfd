## psm4cfd

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)  [![Build Status](https://travis-ci.org/toya42/psm4cfd.svg?branch=master)](https://travis-ci.org/toya42/psm4cfd)  [![](https://github.com/toya42/psm4cfd/workflows/demo/badge.svg)](https://github.com/toya42/psm4cfd/actions)



 psm4cfd is a code of pseudo spectral method  for computational fluid dynamics, wrriten in Fortran.  The psm4cfd can simulate two-dimensional incompressible flow on doubly periodic condition.

<p align="center">
     <img src="https://github.com/toya42/garage/blob/master/psm_01/KH_instability.gif"
width="533" height="480"
alt="Kelvin-Helmholtz instability"
title="Kelvin-Helmholtz instability">
</p>

Kelvin-Helmholtz instability





## Dependency

Fortran2008 (gfortran or ifort)

CMake(3.12 or later)

Intel MKL



## Setup

This section assumes that the OS is Linux. (This program can also be used on macOS.)

### Compiler

#### gfortran

```bash
$sudo apt-get install gfortran
```

#### ifort (Intel® Fortran Compiler)

Intel® Parallel Studio XE

https://software.intel.com/en-us/parallel-studio-xe/choose-download#open-source-contributors



### CMake

* Install

```bash
#install_dir=/install/directory/for/cmake/
$mkdir $install_dir 
$cd $install_dir
$wget https://github.com/Kitware/CMake/releases/download/v3.15.3/cmake-3.15.3.tar.gz  
$tar -xzvf cmake-3.15.3.tar.gz
$mv cmake-3.15.3 cmake_install
$cd cmake_install
$./configure 
$make
$sudo make install
```

Install is completed, add to your "~/.bashrc" file, the following line:

```bash:~/.bashrc
PATH=$install_dir/cmake_install:$install_dir/cmake_install/bin:$PATH
```



### Intel MKL

* Install

  Even if Intel® Parallel Studio XE cannot be used, Intel MKL can be used.

  EURA: https://software.intel.com/en-us/articles/end-user-license-agreement

```bash
$wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
$sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
$sudo wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list 
$sudo apt-get update
$sudo apt-get install intel-mkl-2019.4-070 --allow-unauthenticated
```

Install is completed, add to your "~/.bashrc" file, the following line:

```bash:~/.bashrc
source /opt/intel/mkl/bin/mklvars.sh intel64 lp64
```



## Usage

  According to `demo`, you need to prepare `main.F90` ,`constants.F90` and `initial_flowfield.F90`.

### Compile

 Options of `cmake` are `CMAKE_Fortran_COMPILER`, `CMAKE_BUILD_TYPE` and `large_array`.

You need to type, for example,

```bash
$cmake -D CMAKE_Fortran_COMPILER=ifort    -D CMAKE_BUILD_TYPE=debug -Dlarge_array=OFF
```

or add to your "~/.bashrc" file, the following line:

```bash
# cmake
# small array (DOF<2^14)
#  debug
alias cmakeids='cmake -D CMAKE_Fortran_COMPILER=ifort    -D CMAKE_BUILD_TYPE=debug -Dlarge_array=OFF'
alias cmakegds='cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=debug -Dlarge_array=OFF'
#  fast
alias cmakeifs='cmake -D CMAKE_Fortran_COMPILER=ifort    -D CMAKE_BUILD_TYPE=fast -Dlarge_array=OFF'
alias cmakegfs='cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=fast -Dlarge_array=OFF'

# large array (DOF>2^14)
#  debug
alias cmakeidl='cmake -D CMAKE_Fortran_COMPILER=ifort    -D CMAKE_BUILD_TYPE=debug -Dlarge_array=ON'
alias cmakegdl='cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=debug -Dlarge_array=ON'
#  fast
alias cmakeifl='cmake -D CMAKE_Fortran_COMPILER=ifort    -D CMAKE_BUILD_TYPE=fast -Dlarge_array=ON'
alias cmakegfl='cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_BUILD_TYPE=fast -Dlarge_array=ON'
```

After `$source ~/.bashrc`, you can

```bash
$mkdir build
$cd build
$cmake[???] ..
$make
```

### Execute

```bash
#choose  number of threads (ex. 2)
$export OMP_NUM_THREADS=2
#set stacksize (ex. 1GB)
$export OMP_STACKSIZE=1G
# make output directory
$mkdir output
# execute 
$./2dvorticity_psm.exe
```



### Visulalize

Output file format is plot3d. So any visualization softwares can  be used that supports plot3d format. For example, Paraview, Fieldveiw, etc...

## License

This software is released under the MIT License, see [LICENSE](https://github.com/toya42/psm4cfd/blob/master/LICENSE).

## Authors

toya42

## Contact

Please report bugs and other issues through the issue tracker at:

https://github.com/toya42/psm4cfd/issues

## References

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Fundamentals in Single Domain*, Springer-Verlag (2006).

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Evolution to Complex Geometries and Applications to Fluid Dynamics*, Springer-Verlag (2007).
