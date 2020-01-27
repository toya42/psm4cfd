cmake_minimum_required(VERSION 3.12)

project(linearadvection Fortran)
enable_language(Fortran)

find_package(OpenMP REQUIRED)

set(EXECUTABLE 2dvorticity_psm.exe)


if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set(CMAKE_Fortran_FLAGS_DEBUG   "-O0                -warn all -standf08     \
                                    -warn declarations -CB       -check uninit \
                                    -check all         -fpe0     -traceback    \
                                    -g                 -fpp      -Ddebug -lpthread -ldl")
   set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -fpp -march=native -mtune=native -p -lpthread -ldl")
   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -fpp -march=native -mtune=native  -lpthread -ldl ")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
   set(CMAKE_Fortran_FLAGS_DEBUG   "-O0            -Wall -pedantic  -std=f2008                        \
                                    -fbounds-check -Wuninitialized  -ffpe-trap=invalid,zero,overflow  \
                                    -fbacktrace    -g -Ddebug -lpthread -ldl -Wl,--no-as-needed")
   set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -march=narive -mtune=native -p -lpthread -ldl -Wl,--no-as-needed")
   set(CMAKE_Fortran_FLAGS_FAST    "-O3 -march=native -mtune=native -lpthread -ldl -Wl,--no-as-needed")
endif()

option(large_array "32/64bit integer" OFF)
if(large_array)
   if(APPLE)
      add_definitions(-Dintegertype=1 -i8 -mcmodel=large)
   elseif(UNIX)
      if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
         add_definitions(-Dintegertype=1 -i8 -mcmodel=large)
      elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
         add_definitions(-Dintegertype=1 -fdefault-integer-8 -mcmodel=large)
      endif()
   endif()
else()
   add_definitions(-Dintegertype=0 -mcmodel=large)
endif()

add_executable(${EXECUTABLE}
   main_demo02.F90
   constants_demo02.F90
   initial_flowfield_demo02.F90
   ../../source/fft2d_mkl.F90
   ../../source/arrays_2D_VE_DP.F90
   ../../source/time_integration_2D_VE_DP.F90
   ../../source/output_2D_VE_DP.F90
   ../../source/sort_spectral_coefficients.F90
   ../../source/velocity_2D_VE_DP.F90
   /opt/intel/mkl/include/mkl_dfti.f90
  )

if(APPLE)
   if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
      target_link_libraries(${EXECUTABLE}
         PRIVATE
# Intel(R) Math Kernel Library Link Line Advisor
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
         mkl_intel_thread
         mkl_core
         OpenMP::OpenMP_Fortran
      )
      if(large_array)
         target_link_libraries(${EXECUTABLE}
            PRIVATE
            mkl_intel_ilp64
         )
      else()
         target_link_libraries(${EXECUTABLE}
            PRIVATE
            mkl_intel_lp64
         )
      endif()
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
         message("**** WARNING **********************")
         message(" macOS + gfortran is not supported ")
         message("***********************************")
         message(FATAL_ERROR)
   else()
      message("**** WARNING **********************")
      message(" This compiler is not supported    ")
      message("***********************************")
      message(FATAL_ERROR)
   endif()
elseif(UNIX)
   if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
      target_link_libraries(${EXECUTABLE}
         PRIVATE
# Intel(R) Math Kernel Library Link Line Advisor
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
         mkl_intel_thread
         mkl_core
         OpenMP::OpenMP_Fortran
      )
      if(large_array)
         target_link_libraries(${EXECUTABLE}
            PRIVATE
            mkl_intel_ilp64
         )
      else()
         target_link_libraries(${EXECUTABLE}
            PRIVATE
            mkl_intel_lp64
         )
      endif()
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
      target_link_libraries(${EXECUTABLE}
         PRIVATE
# Intel(R) Math Kernel Library Link Line Advisor
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#        -Wl,--start-group
#        /opt/intel/mkl/lib/intel64/libmkl_gf_lp64.a
#        /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a
#        /opt/intel/mkl/lib/intel64/libmkl_core.a
#        -Wl,--end-group
         /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.so
         /opt/intel/mkl/lib/intel64/libmkl_core.so
         OpenMP::OpenMP_Fortran
      )
      if(large_array)
         target_link_libraries(${EXECUTABLE}
            PRIVATE
            /opt/intel/mkl/lib/intel64/libmkl_gf_ilp64.so
         )
      else()
         target_link_libraries(${EXECUTABLE}
            PRIVATE
            /opt/intel/mkl/lib/intel64/libmkl_gf_lp64.so
         )
      endif()
   else()
      message("**** WARNING **********************")
      message(" This compiler is not supported    ")
      message("***********************************")
      message(FATAL_ERROR)
   endif()
else()
   message("****WARNING**********************")
   message(" This OS is not supported        ")
   message("*********************************")
   message(FATAL_ERROR)
endif()
