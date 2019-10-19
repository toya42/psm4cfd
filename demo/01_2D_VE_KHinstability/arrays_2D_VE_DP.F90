module arrays
   use,intrinsic :: iso_fortran_env
   use constants, only : imax,jmax
   implicit none
   private
   public constant_parameters_initialize
   public cn1ij,cn2ij,cn3ij,cn4ij,cl1ij
   public xi,yj
   public n_length,p_length
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1) :: cn1ij,cn2ij
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1) :: cl1ij,cn3ij,cn4ij
   real(real64),dimension(:),allocatable :: xi,yj
   integer,dimension(2) :: n_length,p_length
   contains
!-----------------------------------------------------------
   subroutine constant_parameters_initialize
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32) :: i,j
         
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               cn2ij(i,j) = -dble(i*i+j*j)
         end do
      end do
      cn2ij(0,0) = 1.0d0
      cn1ij = 1.0d0/cn2ij


      return
   end subroutine constant_parameters_initialize
!-----------------------------------------------------------
end module arrays
