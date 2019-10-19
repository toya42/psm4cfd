module constants
   use,intrinsic :: iso_fortran_env
   implicit none
   real(real64),parameter :: pi = 4.0d0*atan(1d0)
   real(real64),parameter :: nu = 0.0d0

   integer(int32),parameter :: imax = 4*2**6
   integer(int32),parameter :: jmax = imax
   integer(int32),parameter :: ijmax = imax*jmax
   integer(int32),parameter :: ip = imax*3/2
   integer(int32),parameter :: jp = ip
   real(real64),parameter :: dt = 5.0d-3
   real(real64),parameter :: tmax = 20.0d0
!   integer(int32),parameter :: nmax = 4000   ! tmax/dt
!   integer(int32),parameter :: nout = 50
   integer(int32),parameter :: nmax = 10
   integer(int32),parameter :: nout = 10
end module constants

