submodule (flowfield) initial_flowfield_demo01
   use,intrinsic :: iso_fortran_env
!   use constants, only : pi,imax,jmax
   use fft2d_mkl,only : fft2d_execute_forward,fft2d_execute_backward
   implicit none
   contains
!-----------------------------------------------------------
   module subroutine flowfield_initialize_2D_VE_DP(aij,des_n_r2c)
      use mkl_dfti
      use arrays, only : n_length,p_length,
      use sort_spectral_coefficients, only : workc_to_aij
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: aij
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: bij
      type(dfti_descriptor),pointer :: des_n_r2c
      real(real64),dimension(imax,jmax) :: zeta
      real(real64),dimension(imax+2,jmax) :: workc
      integer(int32) :: i,j
      real(real64) :: x,x1,x2,y
      

      allocate(xi(imax))
      allocate(yj(jmax))
      n_length(1) = imax*jmax
      n_length(2) = (imax/2+1)*2*jmax

      do i=1,imax
         xi(i)  = 2.00d0*pi*dble(i-1)/dble(imax)
      end do
      do j=1,jmax
         yj(j)  = 2.00d0*pi*dble(j-1)/dble(jmax)
      end do
      do j=1,jmax
         do i=1,imax
            x  = xi(i)
            y  = yj(j)
            x1 = x-0.50d0*pi
            x2 = x-1.50d0*pi
            zeta(i,j) =4.0d0*(  dexp(-(x1/5.0d-2)**2)*(1+1.0d-2*cos(2.0d0*(y+0.5*pi)))  &
&                              -dexp(-(x2/5.0d-2)**2)*(1+1.0d-2*cos(2.0d0*(y       ))) )
         end do
      end do

      call fft2d_execute_forward(des_n_r2c,n_length,zeta,workc)
      bij = workc_to_aij(imax,jmax,workc)

!$omp parallel workshare
      aij(1,:,:) = bij(1,:,:)*cn1ij(:,:)
      aij(2,:,:) = bij(2,:,:)*cn1ij(:,:)
!$omp end parallel workshare

      return
   end subroutine flowfield_initialize_2D_VE_DP
!-----------------------------------------------------------
end submodule initial_flowfield_demo01

