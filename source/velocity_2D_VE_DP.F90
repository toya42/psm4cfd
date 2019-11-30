module velocity
   use,intrinsic :: iso_fortran_env
   use constants, only : imax,jmax
   use mkl_dfti
   use fft2d_mkl,only : fft2d_execute_forward,fft2d_execute_backward
   implicit none
   private
   public psi_to_velocity!,velocity_to_psi
   contains
!-----------------------------------------------------------
   subroutine psi_to_velocity(des_n_c2r,aij,vel1,vel2)
      use sort_spectral_coefficients, only : aij_to_workc
      use arrays, only : n_length
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      type(dfti_descriptor),pointer :: des_n_c2r
      real(real64),dimension(imax,jmax),intent(out) :: vel1,vel2
#if integertype==0
      integer(int32) :: i,j
#elif integertype==1
      integer(int64) :: i,j
#endif
      real(real64),dimension(imax+2,jmax) :: workc
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: vij
      real(real64) :: dj,di

! u 
!$omp parallel do shared(vij,aij),private(i,j,dj)
      do j=-jmax/2,jmax/2-1
         dj = dble(j)
         do i=0,imax/2-1
            vij(1,i,j) =  dj*aij(2,i,j)
            vij(2,i,j) = -dj*aij(1,i,j)
         end do
      end do
!$omp end parallel do
      workc = aij_to_workc(imax,jmax,vij)
      call fft2d_execute_backward(des_n_c2r,n_length,workc,vel1)

! v
!$omp parallel do shared(vij,aij),private(i,j,di)
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            di = dble(i)
            vij(1,i,j) = -di*aij(2,i,j)
            vij(2,i,j) =  di*aij(1,i,j)
         end do
      end do
!$omp end parallel do
      workc = aij_to_workc(imax,jmax,vij)
      call fft2d_execute_backward(des_n_c2r,n_length,workc,vel2)

   
   end subroutine psi_to_velocity
!-----------------------------------------------------------
end module velocity

