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
      
      real(real64) :: s,as,dkp,dk,ek!,rnd(2)
      complex(real64) :: zeta,iu
      real(real64),dimension(2,0:imax/2,0:jmax/2) :: rnd

      s = 3.0d0
      dkp = 12.0d0
      as = (2.0d0*3.00+1.0d0)**(s+1.0d0)/(8.0d0*6.0d0)
      iu = (0.0d0,1.0d0)

      allocate(xi(imax))
      allocate(yj(jmax))
      do i=1,imax
         xi(i)  = dble(2.00d0*pi*dble(i-1)/dble(imax))
      end do
      do j=1,jmax
         yj(j)  = dble(2.00d0*pi*dble(j-1)/dble(jmax))
      end do

      call random_number(rnd)
      rnd = rnd*2.0d0*pi
      do j=0,jmax/2-1
         do i=0,imax/2-1
            dk = sqrt(dble(i*i+j*j))
            ek = 0.50d0*as/dkp*(dk/dkp)**(2.0d0*s+1.0d0)*exp(-(0.50d0+s)*(dk/dkp)**2)
            zeta = sqrt(dble(dk/pi*ek))*exp(iu*(rnd(1,i,j)+rnd(2,i,j)))
!            bij(1,i,j) = zeta%RE
!            bij(2,i,j) = zeta%IM
            bij(1,i,j) = real(zeta)
            bij(2,i,j) = dimag(zeta)
         end do
      end do
      do j=-jmax/2,-1
         do i=0,imax/2-1
            dk = sqrt(dble(i*i+j*j))
            ek = 0.50d0*as/dkp*(dk/dkp)**(2.0d0*s+1.0d0)*exp(-(0.50d0+s)*(dk/dkp)**2)
            zeta = sqrt(dble(dk/pi*ek))*exp(iu*(rnd(1,i,-j)-rnd(2,i,-j)))
!            bij(1,i,j) = zeta%RE
!            bij(2,i,j) = zeta%IM
            bij(1,i,j) = real(zeta)
            bij(2,i,j) = dimag(zeta)
         end do
      end do

      aij(1,:,:) = bij(1,:,:)*cn1ij(:,:)
      aij(2,:,:) = bij(2,:,:)*cn1ij(:,:)

      return
   end subroutine flowfield_initialize
!-----------------------------------------------------------
end module initial_flowfield

