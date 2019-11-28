module arrays
   use,intrinsic :: iso_fortran_env
   use constants, only : imax,jmax,ip,jp,dt,nu
   implicit none
   private
   public constant_parameters_initialize
   public cn1ij,cn2ij,cn3ij,cn4ij,cl1ij
   public xi,yj
   public n_length,p_length
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1) :: cn1ij,cn2ij
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1) :: cl1ij,cn3ij,cn4ij
   real(real64),dimension(:),allocatable :: xi,yj
#if integertype==0
   integer(int32),dimension(2) :: n_length,p_length
#elif  integertype==1
   integer(int64),dimension(2) :: n_length,p_length
#endif
   contains
!-----------------------------------------------------------
   subroutine constant_parameters_initialize
      use,intrinsic :: iso_fortran_env
      implicit none
#if integertype==0
      integer(int32) :: i,j
#elif  integertype==1
      integer(int64) :: i,j
#endif

! constants for nonlinear term calculation and \zeta<->\psi transformation
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               cn2ij(i,j) = -dble(i*i+j*j)
         end do
      end do
      cn2ij(0,0) = 1.0d0
      cn1ij = 1.0d0/cn2ij

      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               cn3ij(i,j) = dble(i*i-j*j)*cn1ij(i,j)
         end do
      end do
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               cn4ij(i,j) = dble(i*j)*cn1ij(i,j)
         end do
      end do

! constants for linear term calculation
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            cl1ij(i,j) = exp(-0.50d0*dt*nu*dble(i*i+j*j))
         end do
      end do

! integer arrays for fft
      n_length(1) = imax*jmax
      n_length(2) = (imax+2)*jmax
      p_length(1) = ip*jp
      p_length(2) = (ip+2)*jp

      return
   end subroutine constant_parameters_initialize
!-----------------------------------------------------------
end module arrays
