module initial_flowfield
   use,intrinsic :: iso_fortran_env
   use constants, only : pi
   use fft2d_mkl,only : fft2d_execute_forward,fft2d_execute_backward
   implicit none
   contains
!-----------------------------------------------------------
   subroutine flowfield_initialize(imax,jmax,aij,des_n_r2c)
      use mkl_dfti
      use arrays, only : n_length,cn1ij,xi,yj
      use sort_spectral_coefficients, only : workc_to_aij
      implicit none
#if integertype==0
      integer(int32),intent(in) :: imax,jmax
      integer(int32) :: i,j
#elif integertype==1
      integer(int64),intent(in) :: imax,jmax
      integer(int64) :: i,j
#endif
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: aij
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: bij
      type(dfti_descriptor),pointer :: des_n_r2c
      real(real64),dimension(imax,jmax) :: zeta
      real(real64),dimension(imax+2,jmax) :: workc
      real(real64) :: x,x1,x2,y
      
      allocate(xi(imax))
      allocate(yj(jmax))
      do i=1,imax
         xi(i)  = dble(2.00d0*pi*dble(i-1)/dble(imax))
      end do
      do j=1,jmax
         yj(j)  = dble(2.00d0*pi*dble(j-1)/dble(jmax))
      end do
      do j=1,jmax
         do i=1,imax
            x  = xi(i)
            y  = yj(j)
            x1 = x-0.5d0*pi
            x2 = x-1.5d0*pi
            zeta(i,j) =     4.0d0*( exp(-(x1/5.0d-2)**2)*(1.0d0+1.0d-2*cos(2.0d0*(y+0.5d0*pi)))  &
&                                  -exp(-(x2/5.0d-2)**2)*(1.0d0+1.0d-2*cos(2.0d0*(y         ))) ) 
         end do
      end do

      call fft2d_execute_forward(des_n_r2c,n_length,zeta,workc)
      bij = workc_to_aij(imax,jmax,workc)

      aij(1,:,:) = bij(1,:,:)*cn1ij(:,:)
      aij(2,:,:) = bij(2,:,:)*cn1ij(:,:)

      return
   end subroutine flowfield_initialize
!-----------------------------------------------------------
end module initial_flowfield

