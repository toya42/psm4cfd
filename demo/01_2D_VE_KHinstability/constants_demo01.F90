module constants
   use,intrinsic :: iso_fortran_env
   implicit none
!*** user defined constants ***********************************************
! \nu = 1/Re. nu=0 means inviscid calculation
   real(real64),parameter :: nu = 0.0d0
! DOF(physical space & Fourier space). DOF must be a multiple of 4
   integer(int32),parameter :: dof = 4*2**5
! delta t (non-dimensional)
   real(real64),parameter :: dt = 5.0d-3
! loop limit(nmax) & data output interval(nout)
!  --- debug
!   integer(int32),parameter :: nmax = 10
!   integer(int32),parameter :: nout = 10
!  --- demo01
   integer(int32),parameter :: nmax = 4000  
   integer(int32),parameter :: nout = 50
!**************************************************************************

! It is not recommended to change the following constants:
!**************************************************************************
   real(real128),parameter :: pi = 4.0_real128*atan(1.0_real128)
#if integertype==0
   integer(int32),parameter :: imax = dof
   integer(int32),parameter :: jmax = imax
   integer(int32),parameter :: ip = imax*3/2
   integer(int32),parameter :: jp = ip
   integer(int32),parameter :: ijmax = imax*jmax
#elif integertype==1
   integer(int64),parameter :: imax = dof
   integer(int64),parameter :: jmax = imax
   integer(int64),parameter :: ip = imax*3/2
   integer(int64),parameter :: jp = ip
   integer(int64),parameter :: ijmax = imax*jmax
#endif
!**************************************************************************
end module constants

