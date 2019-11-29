## psm4cfd

[![Build Status](https://travis-ci.org/toya42/psm4cfd.svg?branch=master)](https://travis-ci.org/toya42/psm4cfd)  [![](https://github.com/toya42/psm4cfd/workflows/demo/badge.svg)](https://github.com/toya42/psm4cfd/actions)

<p align="center">
     <img src="https://github.com/toya42/garage/blob/master/psm_01/KH_instability.gif"
width="533" height="480"
alt="Kelvin-Helmholtz instability"
title="Kelvin-Helmholtz instability">
</p>

Kelvin-Helmholtz instability





## Dependency

Fortran2008 (gfortran or ifort)

CMake(3.14 or later)

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



## License

This software is released under the MIT License, see [LICENSE](https://github.com/toya42/psm4cfd/blob/master/LICENSE).

## Authors

toya42

Twitter: @toya42_fortran

## References

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Fundamentals in Single Domain*, Springer-Verlag (2006).

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Evolution to Complex Geometries and Applications to Fluid Dynamics*, Springer-Verlag (2007).
