module constants
   use,intrinsic :: iso_fortran_env
   implicit none
!*** user defined constants ***********************************************
! \nu = 1/Re. nu=0 means inviscid calculation
   real(real64),parameter :: nu = 1.0d-3
! DOF(physical space & Fourier space). DOF must be a multiple of 4
   integer(int32),parameter :: dof = 4*2**6
! delta t (non-dimensional)
   real(real64),parameter :: dt = 2.0d-4
! output data (plot3d format)
   integer(int32),parameter :: oflag=6
! 0: no output    
! 1: grid file(.xyz)
! 2: grid file(.xyz) + one   variable(zeta)     function files(.fun)
! 3: grid file(.xyz) + one   variable(zeta)     function files(.fun) + dammy Q file(.q)
! 4: grid file(.xyz) + three variable(zeta,u,v) function files(.fun) 
! 5: grid file(.xyz) + three variable(zeta,u,v) function files(.fun) + dammy Q file(.q)
! 6: grid file(.xyz) + Q files(rho=zeta, u=v=w=e=0) (for paraview)
! loop limit(nmax) & data output interval(nout)
!  --- debug
   integer(int32),parameter :: nmax = 10
   integer(int32),parameter :: nout = 10
!  --- demo01
!   integer(int32),parameter :: nmax = 50000
!   integer(int32),parameter :: nout = 1000
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

