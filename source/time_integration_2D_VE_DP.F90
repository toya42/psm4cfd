! time integration module
!  Governing equation : Two-dimensional vorticity equation
!  Boundary consition : Doubly periodic
!
module time_integration
   use,intrinsic :: iso_fortran_env
   use constants, only : nu,imax,jmax,ip,jp
   use arrays, only : p_length,cl1ij,cn3ij,cn4ij
   implicit none
   private
!   public timeintegration_initialize,rk4_for_nonlinear_ift_for_linear
   public rk4_for_nonlinear_ift_for_linear
   contains
!-----------------------------------------------------------
   subroutine rk4_for_nonlinear_ift_for_linear(dt,ijmax,f0,des_p_r2c,des_p_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      real(real64),intent(in) :: dt
#if integertype==0
      integer(int32),intent(in) :: ijmax
#elif integertype==1
      integer(int64),intent(in) :: ijmax
#endif
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
      real(real64),parameter :: c16 = 1.0d0/6.0d0
      real(real64),parameter :: c13 = 1.0d0/3.0d0
      real(real64),dimension(ijmax),intent(inout) :: f0
      real(real64),dimension(ijmax) :: f1,f2,f3

      call intg_nonlinear(f0,f1,des_p_r2c,des_p_c2r)
!$omp parallel workshare
      f2(:) = f0(:)+0.5d0*dt*f1(:)
      f3(:) = f0(:)+  c16*dt*f1(:)
!$omp end parallel workshare
      call intg_linear(f0)
      call intg_linear(f2)
      call intg_linear(f3)
      call intg_nonlinear(f2,f1,des_p_r2c,des_p_c2r)
      
!$omp parallel workshare
      f2(:) = f0(:)+0.5d0*dt*f1(:)
      f3(:) = f3(:)+  c13*dt*f1(:)
!$omp end parallel workshare
      call intg_nonlinear(f2,f1,des_p_r2c,des_p_c2r)

!$omp parallel workshare
      f2(:) = f0(:)+      dt*f1(:)
      f3(:) = f3(:)+  c13*dt*f1(:)
!$omp end parallel workshare

      call intg_linear(f2)
      call intg_linear(f3)
      call intg_nonlinear(f2,f1,des_p_r2c,des_p_c2r)

!$omp parallel workshare
      f0(:) = f3(:)+  c16*dt*f1(:)
!$omp end parallel workshare

      return
   end subroutine rk4_for_nonlinear_ift_for_linear
!-----------------------------------------------------------
   subroutine intg_nonlinear(f,df,des_p_r2c,des_p_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      use fft2d_mkl, only : fft2d_execute_forward,fft2d_execute_backward
      use sort_spectral_coefficients, only : aij_to_workc,workc_to_aij
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: f
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: df
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
      real(real64),dimension(2,0:ip/2-1,-jp/2:jp/2-1) :: ht_p
      real(real64),dimension(ip+2,jp) :: workc
      real(real64)   ,dimension(0:ip-1,0:jp-1)   :: u,v,uv,v2mu2
#if integertype==0
      integer(int32) :: i,j
#elif integertype==1
      integer(int64) :: i,j
#endif
      real(real64) :: di,dj

! calc. u
!$omp parallel workshare
      ht_p(:,:,:) = 0.0d0
!$omp end parallel workshare
!$omp parallel do shared(ht_p,df),private(i,j,dj)
      do j=-jmax/2,jmax/2-1
         dj = dble(j)
         do i=0,imax/2-1
            ht_p(1,i,j) =  dj*f(2,i,j)
            ht_p(2,i,j) = -dj*f(1,i,j)
         end do
      end do
!$omp end parallel do
      workc = aij_to_workc(ip,jp,ht_p)
      call fft2d_execute_backward(des_p_c2r,p_length,workc,u)

! calc. v
!$omp parallel workshare
      ht_p = 0.0d0
!$omp end parallel workshare
!$omp parallel do shared(ht_p,df),private(j,i,di)
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            di = dble(i)
            ht_p(1,i,j) = -di*f(2,i,j)
            ht_p(2,i,j) =  di*f(1,i,j)
         end do
      end do
!$omp end parallel do
      workc = aij_to_workc(ip,jp,ht_p)
      call fft2d_execute_backward(des_p_c2r,p_length,workc,v)

! calc. Fij
!$omp parallel workshare
      uv(:,:)    = u(:,:)*v(:,:)
!$omp end parallel workshare
      call fft2d_execute_forward(des_p_r2c,p_length,uv,workc)
      ht_p = workc_to_aij(ip,jp,workc)
!$omp parallel do shared(df,cn3ij,ht_p),private(i,j)
       do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            df(1,i,j) = cn3ij(i,j)*ht_p(1,i,j)
            df(2,i,j) = cn3ij(i,j)*ht_p(2,i,j)
         end do
      end do
!$omp end parallel do

!$omp parallel workshare
      v2mu2(:,:) = v(:,:)*v(:,:)-u(:,:)*u(:,:)
!$omp end parallel workshare
      call fft2d_execute_forward(des_p_r2c,p_length,v2mu2,workc)
      ht_p = workc_to_aij(ip,jp,workc)
!$omp parallel do shared(df,cn4ij,ht_p),private(i,j)
       do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            df(1,i,j) = df(1,i,j)+cn4ij(i,j)*ht_p(1,i,j)
            df(2,i,j) = df(2,i,j)+cn4ij(i,j)*ht_p(2,i,j)
         end do
      end do
!$omp end parallel do
    
      return

   end subroutine intg_nonlinear
!-----------------------------------------------------------
   subroutine intg_linear(f)
      use,intrinsic :: iso_fortran_env
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(inout) :: f
#if integertype==0
      integer(int32) :: i,j
#elif integertype==1
      integer(int64) :: i,j
#endif

!$omp parallel do shared(f,cl1ij),private(i,j)
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            f(1,i,j) = f(1,i,j)*cl1ij(i,j)
            f(2,i,j) = f(2,i,j)*cl1ij(i,j)
         end do
      end do
!$omp end parallel do
      return

   end subroutine intg_linear
!-----------------------------------------------------------
end module time_integration

