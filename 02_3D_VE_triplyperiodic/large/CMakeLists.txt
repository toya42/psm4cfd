cmake_minimum_required(VERSION 3.14)

project(linearadvection Fortran)
enable_language(Fortran)

find_package(OpenMP REQUIRED)

set(EXECUTABLE 3dvorticity_psm.exe)


if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set(CMAKE_Fortran_FLAGS_DEBUG   "-O0                -warn all -standf08     \
                                    -warn declarations -CB       -check uninit \
                                    -check all         -fpe0     -traceback    \
                                    -g                 -fpp      -Ddebug       \
                                    -check bounds -mcmodel=large -i8")
   set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -fpp -march=native -mtune=native -p")
#  set(CMAKE_Fortran_FLAGS_FAST    "-O3 -fpp -march=native -mtune=native -mkl -static-intel -liomp5 -lpthread -ldl")
#   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -fpp -march=native -mtune=native -lpthread -ldl -mcmodel=large -Ddebug -i8")
   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -fpp -march=native -mtune=native -lpthread -ldl -mcmodel=large -i8")
#   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -fpp -march=native -mtune=native -mcmodel=large")
#   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -fpp   -mcmodel=large")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
   set(CMAKE_Fortran_FLAGS_DEBUG   "-O0            -Wall -pedantic  -std=f2008                        \
                                    -fbounds-check -Wuninitialized  -ffpe-trap=invalid,zero,overflow  \
                                    -fbacktrace    -g -Ddebug")
   set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -march=narive -mtune=native -p")
   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -march=native -mtune=native")
endif()

add_executable(${EXECUTABLE}
   main.F90
   /opt/intel/mkl/include/mkl_dfti.f90
  )

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   target_link_libraries(${EXECUTABLE}
      PRIVATE
# Intel Link Advisor
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#=== A. static library
      -Wl,--start-group
#      /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a
      /opt/intel/mkl/lib/intel64/libmkl_intel_ilp64.a
      /opt/intel/mkl/lib/intel64/libmkl_intel_thread.a
#      /opt/intel/mkl/lib/intel64/libmkl_sequential.a
      /opt/intel/mkl/lib/intel64/libmkl_core.a
      -Wl,--end-group
##=== B. shared library
##      mkl_intel_lp64
#      mkl_intel_ilp64
##      mkl_intel_thread
#      mkl_sequential
#      mkl_core

      OpenMP::OpenMP_Fortran
   )
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
   target_link_libraries(${EXECUTABLE}
      PRIVATE
# Intel Link Advisor
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
      mkl_gf_lp64
      mkl_gnu_thread
      mkl_core
      OpenMP::OpenMP_Fortran
   )
endif()


