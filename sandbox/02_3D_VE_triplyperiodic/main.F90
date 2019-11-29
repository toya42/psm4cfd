!===========================================================
! governing equation : 3d vorticity equation
! bondary condition : triply periodic
! numerical method
!     spatial discretization : spectral method
!        trial and test function : Fourier series
!     non-linear term : pseudo spectral method
!        de-aliasing : 3/2-rule (padding)
!     time integration : Runge-Kutta 4 (nonlinear term) 
!                        integrating factor (linear term)
!===========================================================
!                                2019.09.24 : ver.1
!===========================================================
module constants
   use,intrinsic :: iso_fortran_env
   implicit none
   real(real64),parameter :: pi = 4.0d0*atan(1d0)
!   real(real64),parameter :: nu = 1.0d-2     ! 1/Re
   real(real64),parameter :: nu = 6.25d-4     ! 1/Re
   real(real64),parameter :: gmmai = 1.0d0/1.40d0

   integer(int32),parameter :: imax = 4*2**5
   integer(int32),parameter :: jmax = imax
   integer(int32),parameter :: kmax = imax
   integer(int32),parameter :: ip = imax*3/2
   integer(int32),parameter :: jp = ip
   integer(int32),parameter :: kp = ip
   real(real64),parameter :: dt = 1.0d-3
   real(real64),parameter :: tmax = 10.0d0
   integer(int32),parameter :: nmax = 8000   ! tmax/dt
   integer(int32),parameter :: nout = 250
!   integer(int32),parameter :: nmax = 1
!   integer(int32),parameter :: nout = 1
end module constants
!===========================================================
module arrays
   use,intrinsic :: iso_fortran_env
   use constants
   implicit none
   real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3) :: voijk
end module arrays
!===========================================================
module fft_mkl
   use,intrinsic :: iso_fortran_env
   use mkl_dfti
   implicit none
   private
   public fft_initialize,fft_finalize,fft_execute_forward,fft_execute_backward
   contains
!-----------------------------------------------------------
   subroutine fft_initialize(imax,jmax,kmax,des_r2c,des_c2r)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax,kmax
      type(dfti_descriptor),pointer :: des_r2c,des_c2r
      integer(int32),dimension(12) :: dftistatus
      integer(int32),dimension(3) :: length
      integer(int32),dimension(4) :: stride_i,stride_o
      integer(int32) :: i
      real(real64) :: ni

      ni = 1.0d0/dble(imax*jmax*kmax)
      length(1) = imax
      length(2) = jmax
      length(3) = kmax
      stride_i(1) = 0
      stride_i(2) = 1
      stride_i(3) = imax
      stride_i(4) = imax*jmax
      stride_o(1) = 0
      stride_o(2) = 1
      stride_o(3) = imax/2+1
      stride_o(4) = (imax/2+1)*jmax
!r2c
      dftistatus = 0
      dftistatus( 1) = dfticreatedescriptor(des_r2c,dfti_double,dfti_real,3,length)
      dftistatus( 2) = dftisetvalue        (des_r2c,dfti_forward_scale,ni)
      dftistatus( 3) = dftisetvalue        (des_r2c,dfti_placement,dfti_not_inplace)
      dftistatus( 4) = dftisetvalue        (des_r2c,dfti_conjugate_even_storage,dfti_complex_complex)
      dftistatus( 5) = dftisetvalue        (des_r2c,dfti_packed_format,dfti_cce_format)
      dftistatus( 6) = dftisetvalue        (des_r2c,dfti_input_strides, stride_i)
      dftistatus( 7) = dftisetvalue        (des_r2c,dfti_output_strides,stride_o)
      dftistatus( 8) = dfticommitdescriptor(des_r2c)
!c2r
      dftistatus( 9) = dfticopydescriptor  (des_r2c,des_c2r)
      dftistatus(10) = dftisetvalue        (des_c2r,dfti_input_strides, stride_o)
      dftistatus(11) = dftisetvalue        (des_c2r,dfti_output_strides,stride_i)
      dftistatus(12) = dfticommitdescriptor(des_c2r)

      print *,'fft_initialize) imax   =',imax
      print *,'fft_initialize) jmax   =',jmax
      print *,'fft_initialize) kmax   =',kmax
#ifdef debug
      do i=1,12
         if(dftistatus(i).ne.0) then
            print *,'fft_initialize) dft setting error:',i,dftistatus(i)
            stop
         end if
      end do
      print *,'fft_initialize) dft settings are completed'
#endif

      return
   end subroutine fft_initialize
!-----------------------------------------------------------
   subroutine fft_finalize(des_r2c,des_c2r)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_r2c,des_c2r
      integer(int32),dimension(2) :: dftistatus
      integer(int32) :: i

      dftistatus(1)=DftiFreeDescriptor(Des_r2c)
      dftistatus(2)=DftiFreeDescriptor(Des_c2r)

#ifdef debug
      do i=1,2
         if(dftistatus(i).ne.0) then
            print *,'fft_finalize) dft finalization error:',i
            stop
         end if
      end do
      print *,'fft_finalize) dft finalization is completed'
#endif
      return
   end subroutine fft_finalize
!-----------------------------------------------------------
   subroutine fft_execute_forward(des_r2c,length,real3d,complex3d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_r2c
      integer(int32),dimension(2),intent(in) :: length
      real(real64),dimension(length(1)),intent(in)  :: real3d
      real(real64),dimension(length(2)),intent(out) :: complex3d
      integer(int32) :: dftistatus

      dftistatus = dfticomputeforward(des_r2c,real3d,complex3d)
#ifdef debug
      if(dftistatus.ne.0) then
         print *,'fft_execute_forward) dft finalization error:',dftistatus
         stop
      end if
#endif
      return
   end subroutine fft_execute_forward
!-----------------------------------------------------------
   subroutine fft_execute_backward(des_c2r,length,complex3d,real3d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_c2r
      integer(int32),dimension(2),intent(in) :: length
      real(real64),dimension(length(2)),intent(in)  :: complex3d
      real(real64),dimension(length(1)),intent(out) :: real3d
      integer(int32) :: dftistatus

      dftistatus = dfticomputebackward(des_c2r,complex3d,real3d)
#ifdef debug
      if(dftistatus.ne.0) then
         print *,'fft_execute_forward) dft finalization error:',dftistatus
         stop
      end if
#endif
      return
   end subroutine fft_execute_backward
!-----------------------------------------------------------
end module fft_mkl
!===========================================================
module flowfield
   use,intrinsic :: iso_fortran_env
   use constants, only : pi,imax,jmax,kmax,ip,jp,kp
   use fft_mkl,only : fft_execute_forward,fft_execute_backward
   implicit none
   private
   public flowfield_initialize,wrtd
   public aijk_to_workc,workc_to_aijk
   public velocity_to_vorticity,vorticity_to_velocity,velocity_to_pressure
   public guarantee_solenoidal
   public cn1ijk,cn2ijk
   public n_length,p_length
   real(real64),dimension(:),allocatable :: xi,yj,zk
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1) :: cn1ijk,cn2ijk
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,6) :: cn3ijk
   integer(int32),dimension(2) :: n_length,p_length
   contains
!-----------------------------------------------------------
   subroutine flowfield_initialize(voijk,des_n_r2c)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(out) :: voijk
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3) :: veijk
      type(dfti_descriptor),pointer :: des_n_r2c
      real(real64),dimension(imax,jmax,kmax) :: u,v,w
      real(real64),dimension(imax+2,jmax,kmax) :: workc
      integer(int32) :: i,j,k
      real(real64) :: x,y,z
      
! constants for aij<->bij and  nonlinear term calculation
!  #1 and #2
      do k = -kmax/2,kmax/2-1
         do j = -jmax/2,jmax/2-1
            do i = 0,imax/2-1
               cn2ijk(i,j,k) = dble(i*i+j*j+k*k)
            end do
         end do
      end do
      cn2ijk(0,0,0) = 1.0d0
      cn1ijk = 1.0d0/cn2ijk
      cn2ijk(0,0,0) = 0.0d0

      do k = -kmax/2,kmax/2-1
         do j = -jmax/2,jmax/2-1
            do i = 0,imax/2-1
               cn3ijk(i,j,k,1) = -cn1ijk(i,j,k)*dble(i*i)
               cn3ijk(i,j,k,2) = -cn1ijk(i,j,k)*dble(j*j)
               cn3ijk(i,j,k,3) = -cn1ijk(i,j,k)*dble(k*k)
               cn3ijk(i,j,k,4) = -cn1ijk(i,j,k)*dble(i*j)
               cn3ijk(i,j,k,5) = -cn1ijk(i,j,k)*dble(j*k)
               cn3ijk(i,j,k,6) = -cn1ijk(i,j,k)*dble(k*i)
            end do
         end do
      end do


      allocate(xi(imax))
      allocate(yj(jmax))
      allocate(zk(kmax))
      n_length(1) = imax*jmax*kmax
      n_length(2) = (imax+2)*jmax*kmax
      p_length(1) = ip*jp*kp
      p_length(2) = (ip+2)*jp*kp   ! (ip+2)*jp

      do i = 1,imax
         xi(i) = 2.00d0*pi*dble(i-1)/dble(imax)
      end do
      do j = 1,jmax
         yj(j) = 2.00d0*pi*dble(j-1)/dble(jmax)
      end do
      do k = 1,kmax
         zk(k) = 2.00d0*pi*dble(k-1)/dble(kmax)
      end do

! Taylor-Green vortex
      do k = 1,kmax
         do j = 1,jmax
            do i = 1,imax
               x = xi(i)
               y = yj(j)
               z = zk(k)
! definition of International workshop on high order CFD methods
               u(i,j,k) = sin(x-pi)*cos(y-pi)*cos(z-pi)
               v(i,j,k) =-cos(x-pi)*sin(y-pi)*cos(z-pi)
               w(i,j,k) = 0.0d0
            end do
         end do
      end do


! velocity: pysical(u,v,w) to fourier(veijk)
!-- u
      call fft_execute_forward(des_n_r2c,n_length,u,workc)
      veijk(:,:,:,:,1) = workc_to_aijk(imax,jmax,kmax,workc)
!-- v
      call fft_execute_forward(des_n_r2c,n_length,v,workc)
      veijk(:,:,:,:,2) = workc_to_aijk(imax,jmax,kmax,workc)
!-- w
      call fft_execute_forward(des_n_r2c,n_length,w,workc)
      veijk(:,:,:,:,3) = workc_to_aijk(imax,jmax,kmax,workc)

! fourier coefficient of vorticity (aij)
      call velocity_to_vorticity(veijk,voijk)

! Are vorticity and velocity solenoidal?
      call checksolenoidal(veijk,voijk,-1)

      return
   end subroutine flowfield_initialize
!-----------------------------------------------------------
   subroutine velocity_to_vorticity(veijk,voijk)
      use,intrinsic :: iso_fortran_env
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(in)  :: veijk
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(out) :: voijk
      integer(int32) :: i,j,k
      real(real64) :: di,dj,dk
      real(real64),dimension(2) :: uht,vht,wht

!$omp parallel do shared(veijk,voijk),private(i,j,k,di,dj,dk,uht,vht,wht)
      do k = -kmax/2,kmax/2-1
         dk  = dble(k)
         do j = -jmax/2,jmax/2-1
            dj  = dble(j)
            do i = 0,imax/2-1
               di =dble(i)
               uht(1) = veijk(1,i,j,k,1)
               uht(2) = veijk(2,i,j,k,1)
               vht(1) = veijk(1,i,j,k,2)
               vht(2) = veijk(2,i,j,k,2)
               wht(1) = veijk(1,i,j,k,3)
               wht(2) = veijk(2,i,j,k,3)
! w_1=Dw/Dy-Dv/Dz
               voijk(1,i,j,k,1) = -dj*wht(2)+dk*vht(2)
               voijk(2,i,j,k,1) =  dj*wht(1)-dk*vht(1)
! w_2=Du/Dz-Dw/Dx
               voijk(1,i,j,k,2) = -dk*uht(2)+di*wht(2)
               voijk(2,i,j,k,2) =  dk*uht(1)-di*wht(1)
! w_3=Dv/Dx-Du/Dy
               voijk(1,i,j,k,3) = -di*vht(2)+dj*uht(2)
               voijk(2,i,j,k,3) =  di*vht(1)-dj*uht(1)
            end do
         end do
      end do
!$omp end parallel do


      return
   end subroutine velocity_to_vorticity
!-----------------------------------------------------------
   subroutine vorticity_to_velocity(voijk,veijk)
      use,intrinsic :: iso_fortran_env
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(in)  :: voijk
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(out) :: veijk
      integer(int32) :: i,j,k
      real(real64) :: di,dj,dk
      real(real64),dimension(2) :: o1ht,o2ht,o3ht  ! \hat{Omega}
      real(real64) :: ck

!$omp parallel do shared(voijk,veijk,cn1ijk),private(i,j,k,di,dj,dk,o1ht,o2ht,o3ht,ck)
      do k = -kmax/2,kmax/2-1
         dk  = dble(k)
         do j = -jmax/2,jmax/2-1
            dj  = dble(j)
            do i = 0,imax/2-1
               di = dble(i)
               ck = cn1ijk(i,j,k)
               o1ht(1) = voijk(1,i,j,k,1)
               o1ht(2) = voijk(2,i,j,k,1)
               o2ht(1) = voijk(1,i,j,k,2)
               o2ht(2) = voijk(2,i,j,k,2)
               o3ht(1) = voijk(1,i,j,k,3)
               o3ht(2) = voijk(2,i,j,k,3)
! uht
               veijk(1,i,j,k,1) = (-dj*o3ht(2)+dk*o2ht(2))*ck
               veijk(2,i,j,k,1) = ( dj*o3ht(1)-dk*o2ht(1))*ck
! vht
               veijk(1,i,j,k,2) = (-dk*o1ht(2)+di*o3ht(2))*ck
               veijk(2,i,j,k,2) = ( dk*o1ht(1)-di*o3ht(1))*ck
! wht
               veijk(1,i,j,k,3) = (-di*o2ht(2)+dj*o1ht(2))*ck
               veijk(2,i,j,k,3) = ( di*o2ht(1)-dj*o1ht(1))*ck
            end do
         end do
      end do
!$omp end parallel do

      return
   end subroutine vorticity_to_velocity
!-----------------------------------------------------------
   subroutine velocity_to_pressure(veijk,pijk,des_p_r2c,des_p_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      use fft_mkl, only : fft_execute_forward,fft_execute_backward
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(in)  :: veijk
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1  ),intent(out) :: pijk
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
      real(real64),dimension(ip+2,jp,kp) :: workc
      real(real64),dimension(2,0:ip/2-1,-jp/2:jp/2-1,-kp/2:kp/2-1) :: ht_p
      real(real64),dimension(0:ip-1,0:jp-1,0:kp-1) :: u,v,w
      real(real64),dimension(0:ip-1,0:jp-1,0:kp-1,6) :: uu
      integer(int32) :: i,j,k,n

!$omp parallel workshare
      pijk(:,:,:,:) = 0.0d0
!$omp end parallel workshare


! \hat{velocity} to velocity  using  IFFT with 3/2 rule
! u
!$omp parallel workshare
      ht_p(:,:,:,:) = 0.0d0
!$omp end parallel workshare
!$omp parallel workshare
      ht_p(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1)  &
&          =veijk(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,1)
!$omp end parallel workshare
      workc = aijk_to_workc(ip,jp,kp,ht_p)
      call fft_execute_backward(des_p_c2r,p_length,workc,u)
! v
!$omp parallel workshare
      ht_p(:,:,:,:) = 0.0d0
!$omp end parallel workshare
!$omp parallel workshare
      ht_p(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1)  &
&          =veijk(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,2)
!$omp end parallel workshare
      workc = aijk_to_workc(ip,jp,kp,ht_p)
      call fft_execute_backward(des_p_c2r,p_length,workc,v)
! v
!$omp parallel workshare
      ht_p(:,:,:,:) = 0.0d0
!$omp end parallel workshare
!$omp parallel workshare
      ht_p(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1)  &
&          =veijk(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3)
!$omp end parallel workshare
      workc = aijk_to_workc(ip,jp,kp,ht_p)
      call fft_execute_backward(des_p_c2r,p_length,workc,w)

! calc. u_i*u_j
!$omp parallel workshare
      uu(:,:,:,1) = u(:,:,:)*u(:,:,:)
      uu(:,:,:,2) = v(:,:,:)*v(:,:,:)
      uu(:,:,:,3) = w(:,:,:)*w(:,:,:)
      uu(:,:,:,4) = u(:,:,:)*v(:,:,:)*2.0d0
      uu(:,:,:,5) = v(:,:,:)*w(:,:,:)*2.0d0
      uu(:,:,:,6) = w(:,:,:)*u(:,:,:)*2.0d0
!$omp end parallel workshare

! u_i*u_j to \widetilde{u_i*u_j} using FFT with 3/2 rule -> pijk (\hat{p})
      do n=1,6
         call fft_execute_forward(des_p_r2c,p_length,uu(:,:,:,n),workc)
         ht_p = workc_to_aijk(ip,jp,kp,workc)
!$omp parallel do shared(pijk,ht_p,cn3ijk),private(i,j,k)
         do k = -kmax/2,kmax/2-1
            do j = -jmax/2,jmax/2-1
               do i = 0,imax/2-1
                  pijk(1,i,j,k)=pijk(1,i,j,k)+ht_p(1,i,j,k)*cn3ijk(i,j,k,n)
                  pijk(2,i,j,k)=pijk(2,i,j,k)+ht_p(2,i,j,k)*cn3ijk(i,j,k,n)
               end do
            end do
         end do
!$omp end parallel do
      end do


      return
   end subroutine velocity_to_pressure
!-----------------------------------------------------------
   subroutine guarantee_solenoidal(voijk)
      use,intrinsic :: iso_fortran_env
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(inout) :: voijk ! vorticity
      integer(int32) :: i,j,k
      real(real64) :: di,dj,dk,dki,dji

!$omp parallel do shared(voijk),private(i,j,k,di,dj,dk,dki)
      do k = -kmax/2,-1
         dk  = dble(k)
         dki = 1.0d0/dk
         do j = -jmax/2,jmax/2-1
            dj  = dble(j)
            do i = 0,imax/2-1
               di = dble(i)
               voijk(1,i,j,k,3) = -dki*(di*voijk(1,i,j,k,1)+dj*voijk(1,i,j,k,2))
               voijk(2,i,j,k,3) = -dki*(di*voijk(2,i,j,k,1)+dj*voijk(2,i,j,k,2))
            end do
         end do
      end do
!$omp end parallel do

      k=0
!$omp parallel do shared(voijk),private(i,j,di,dj,dji)
      do j = -jmax/2,-1
         dj  = dble(j)
         dji = 1.0d0/dj
         do i = 0,imax/2-1
            di = dble(i)
            voijk(1,i,j,k,2) = -dji*di*voijk(1,i,j,k,1)
            voijk(2,i,j,k,2) = -dji*di*voijk(2,i,j,k,1)
         end do
      end do
!$omp end parallel do
      j=0
      i=0
      voijk(1,i,j,k,3) = 0.0d0
      voijk(2,i,j,k,3) = 0.0d0
!$omp parallel do shared(voijk),private(i,j,di,dj,dji)
      do j = 1,jmax/2-1
         dj  = dble(j)
         dji = 1.0d0/dj
         do i = 0,imax/2-1
            di = dble(i)
            voijk(1,i,j,k,2) = -dji*di*voijk(1,i,j,k,1)
            voijk(2,i,j,k,2) = -dji*di*voijk(2,i,j,k,1)
         end do
      end do
!$omp end parallel do    
!$omp parallel do shared(voijk),private(i,j,k,di,dj,dk,dki)
      do k = 1,kmax/2-1
         dk  = dble(k)
         dki = 1.0d0/dk
         do j = -jmax/2,jmax/2-1
            dj  = dble(j)
            do i = 0,imax/2-1
               di = dble(i)
               voijk(1,i,j,k,3) = -dki*(di*voijk(1,i,j,k,1)+dj*voijk(1,i,j,k,2))
               voijk(2,i,j,k,3) = -dki*(di*voijk(2,i,j,k,1)+dj*voijk(2,i,j,k,2))
            end do
         end do
      end do
!$omp end parallel do

     return
   end subroutine guarantee_solenoidal

!-----------------------------------------------------------
   subroutine checksolenoidal(veijk,voijk,loop)
      use,intrinsic :: iso_fortran_env
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(in) :: veijk ! velocity
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(in) :: voijk ! vorticity
      integer(int32),intent(in) :: loop
      real(real64) :: re,im,di,dj,dk
      integer(int32) :: i,j,k

!$omp parallel do shared(veijk,voijk,loop),private(i,j,k,di,dj,dk,re,im)
       do k = -kmax/2,kmax/2-1
         dk = dble(k)
         do j = -jmax/2,jmax/2-1
            dj = dble(j)
            do i = 0,imax/2-1
               di = dble(i)
               re = di*veijk(1,i,j,k,1)+dj*veijk(1,i,j,k,2)+dk*veijk(1,i,j,k,3)
               im = di*veijk(2,i,j,k,1)+dj*veijk(2,i,j,k,2)+dk*veijk(2,i,j,k,3)
               if(abs(re)>1.0d-15 .or. abs(im)>1.0d-15) then
                  write(20,'(4i8,2E18.8e3)'),loop,i,j,k,re,im
                  print *,'checksolenoidal) velocity error:',i,j,k
               end if
               re = di*voijk(1,i,j,k,1)+dj*voijk(1,i,j,k,2)+dk*voijk(1,i,j,k,3)
               im = di*voijk(2,i,j,k,1)+dj*voijk(2,i,j,k,2)+dk*voijk(2,i,j,k,3)
               if(abs(re)>1.0d-15 .or. abs(im)>1.0d-15) then
                  write(21,'(4i8,2E18.8e3)'),loop,i,j,k,re,im
                  print *,'checksolenoidal) vorticity error:',i,j,k
               end if
            end do
         end do
      end do
!$omp end parallel do

      return
   end subroutine checksolenoidal
!-----------------------------------------------------------
   subroutine wrtd(loop,voijk,des_n_c2r,des_p_r2c,des_p_c2r)
      use,intrinsic :: iso_fortran_env
      use constants, only : gmmai,nu,dt
      use mkl_dfti
      use fft_mkl, only : fft_execute_backward
      implicit none
      integer(int32),intent(in) :: loop
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(in) :: voijk ! vorticity
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3) :: veijk            ! velocity
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1)   :: pijk             ! pressure
      type(dfti_descriptor),pointer :: des_n_c2r,des_p_r2c,des_p_c2r
      real(real64),dimension(imax,jmax,kmax,3) :: velocity,vorticity
      real(real64),dimension(imax,jmax,kmax) :: pressure,q_criterion

      integer(int32) :: i,j,k,n
      real(real64),dimension(imax+2,jmax,kmax) :: workc
      real(real64),dimension(:,:,:,:),allocatable :: buf
      character(len=50) :: filename
      real(real64) :: dammy,u,v,w,q,ek,dekdt,enstrophy

! output grid file (./output/grid.xyz) 
      if(loop.eq.0) then
         allocate(buf(imax,jmax,kmax,3))
         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  buf(i,j,k,1) = xi(i)
                  buf(i,j,k,2) = yj(j)
                  buf(i,j,k,3) = zk(k)
               end do
            end do
         end do
         open(10,file='output/grid.xyz',form='unformatted',access='stream',status='replace')
         i=imax
         j=jmax
         k=kmax
         write(10) i,j,k
         write(10) buf
         close(10)
         deallocate(buf)
         deallocate(xi,yj,zk)
      end if

! vorticity
      do n=1,3
         workc = aijk_to_workc(imax,jmax,kmax,voijk(:,:,:,:,n))
         call fft_execute_backward(des_n_c2r,n_length,workc,vorticity(:,:,:,n))
      end do

! vorticity to velocity
      call vorticity_to_velocity(voijk,veijk)
! velocity
      do n=1,3
         workc = aijk_to_workc(imax,jmax,kmax,veijk(:,:,:,:,n))
         call fft_execute_backward(des_n_c2r,n_length,workc,velocity(:,:,:,n))
      end do

! pressure
      call velocity_to_pressure(veijk,pijk,des_p_r2c,des_p_c2r)
      workc = aijk_to_workc(imax,jmax,kmax,pijk)
      call fft_execute_backward(des_n_c2r,n_length,workc,pressure)
! Q-criterion
!$omp parallel workshare
      pijk(1,:,:,:) = -pijk(1,:,:,:)*cn2ijk(:,:,:)
      pijk(2,:,:,:) = -pijk(2,:,:,:)*cn2ijk(:,:,:)
!$omp end parallel workshare
      workc = aijk_to_workc(imax,jmax,kmax,pijk)
      call fft_execute_backward(des_n_c2r,n_length,workc,q_criterion)
!$omp parallel workshare
      q_criterion(:,:,:)=0.50d0*q_criterion(:,:,:)
!$omp end parallel workshare

! error estimation (velocity & vorticity must be solenoidal)
!      call checksolenoidal(veijk,voijk,loop)


! q file
      ek = 0.0d0
      enstrophy = 0.0d0
      allocate(buf(imax,jmax,kmax,5))
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
!               u=velocity(i,j,k,1)
!               v=velocity(i,j,k,2)
!               w=velocity(i,j,k,3)
               u=vorticity(i,j,k,1)
               v=vorticity(i,j,k,2)
               w=vorticity(i,j,k,3)
               q=u*u+v*v+w*w
               enstrophy = enstrophy + q
               q=sqrt(q)
               u=velocity(i,j,k,1)
               v=velocity(i,j,k,2)
               w=velocity(i,j,k,3)
               ek = ek + u*u+v*v+w*w
               buf(i,j,k,1)=q
               buf(i,j,k,2)=u*q
               buf(i,j,k,3)=v*q
               buf(i,j,k,4)=w*q
!               buf(i,j,k,1)=1.0d0
!               buf(i,j,k,2)=u
!               buf(i,j,k,3)=v
!               buf(i,j,k,4)=w
               buf(i,j,k,5)=pressure(i,j,k)*gmmai+0.50d0*(u*u+v*v+w*w)
            end do
         end do
      end do
      write(filename,'("output/tgv",i5.5,".q")') loop
      open(10,file=filename,form='unformatted',access='stream',status='replace')
      i=imax
      j=jmax
      k=kmax
      dammy = 0.0d0
      write(10) i,j,k
      write(10) dammy,dammy,dammy,dammy
      write(10) buf
      close(10)
      deallocate(buf)

      ek = 0.50d0*ek/dble(imax*jmax*kmax)
      enstrophy = 0.50d0*enstrophy/dble(imax*jmax*kmax)
      dekdt = 2.0d0*nu*enstrophy
      write(22,'(4E18.8e3)'),dt*dble(loop),ek,dekdt,enstrophy


! function file
      allocate(buf(imax,jmax,kmax,5))
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               buf(i,j,k,1)=vorticity(i,j,k,1)
               buf(i,j,k,2)=vorticity(i,j,k,2)
               buf(i,j,k,3)=vorticity(i,j,k,3)
               buf(i,j,k,4)=pressure(i,j,k)
               buf(i,j,k,5)=q_criterion(i,j,k)
            end do
         end do
      end do
      write(filename,'("output/tgv",i5.5,".fun")') loop
      open(10,file=filename,form='unformatted',access='stream',status='replace')
      i=imax
      j=jmax
      k=kmax
      n=5
      write(10) i,j,k,n
      write(10) buf
      close(10)
      deallocate(buf)


      print *,'wrtd) flowfield output,loop=',loop

      return
   end subroutine wrtd
!-----------------------------------------------------------
   function workc_to_aijk(imax,jmax,kmax,workc) result(aijk)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax,kmax
      real(real64),dimension(imax+2,jmax,kmax),intent(in) :: workc
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1) :: aijk
      integer(int32) :: i,i2,j,j1,jmjm,k,k1,kmkm

!$omp parallel do shared(aijk,workc),private(i,i2,j,j1,jmjm,k,k1)
      do k=0,kmax/2-1
         k1 = k+1
         do j=0,jmax/2-1
            j1 = j+1
            do i=0,imax/2-1
               i2 = i+i
               aijk(1,i,j,k) = workc(i2+1,j1,k1)
               aijk(2,i,j,k) = workc(i2+2,j1,k1)
            end do
         end do
         do j=jmax/2,jmax-1
            j1 = j+1
            jmjm = j-jmax
            do i=0,imax/2-1
               i2 = i+i
               aijk(1,i,jmjm,k) = workc(i2+1,j1,k1) 
               aijk(2,i,jmjm,k) = workc(i2+2,j1,k1) 
            end do
         end do
      end do
!$omp end parallel do
!$omp parallel do shared(aijk,workc),private(i,i2,j,j1,jmjm,k,k1,kmkm)
      do k=kmax/2,kmax-1
         k1 = k+1
         kmkm = k-kmax
         do j=0,jmax/2-1
            j1 = j+1
            do i=0,imax/2-1
               i2 = i+i
               aijk(1,i,j,kmkm) = workc(i2+1,j1,k1)
               aijk(2,i,j,kmkm) = workc(i2+2,j1,k1)
            end do
         end do
         do j=jmax/2,jmax-1
            j1 = j+1
            jmjm = j-jmax
            do i=0,imax/2-1
               i2 = i+i
               aijk(1,i,jmjm,kmkm) = workc(i2+1,j1,k1) 
               aijk(2,i,jmjm,kmkm) = workc(i2+2,j1,k1) 
            end do
         end do
      end do
!$omp end parallel do
      
   end function workc_to_aijk
!-----------------------------------------------------------
   function aijk_to_workc(imax,jmax,kmax,aijk) result(workc)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax,kmax
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1),intent(in) :: aijk
      real(real64),dimension(imax+2,jmax,kmax) :: workc
      integer(int32) :: i,i2,j,j1,jmjm,k,k1,kmkm

!$omp parallel do shared(aijk,workc),private(i,i2,j,j1,jmjm,k,k1)
      do k=0,kmax/2-1
         k1 = k+1
         do j=0,jmax/2-1
            j1 = j+1
            do i=0,imax/2-1
               i2 = i+i
               workc(i2+1,j1,k1) = aijk(1,i,j,k)
               workc(i2+2,j1,k1) = aijk(2,i,j,k)
            end do
            workc(imax+1,j1,k1) = 0.0d0
            workc(imax+2,j1,k1) = 0.0d0
         end do
         do j=jmax/2,jmax-1
            j1 = j+1
            jmjm = j-jmax
            do i=0,imax/2-1
               i2 = i+i
               workc(i2+1,j1,k1) = aijk(1,i,jmjm,k)
               workc(i2+2,j1,k1) = aijk(2,i,jmjm,k)
            end do
            workc(imax+1,j1,k1) = 0.0d0
            workc(imax+2,j1,k1) = 0.0d0
         end do
      end do
!$omp end parallel do
!$omp parallel do shared(aijk,workc),private(i,i2,j,j1,jmjm,k,k1,kmkm)
      do k=kmax/2,kmax-1
         k1 = k+1
         kmkm = k-kmax
         do j=0,jmax/2-1
            j1 = j+1
            do i=0,imax/2-1
               i2 = i+i
               workc(i2+1,j1,k1) = aijk(1,i,j,kmkm)
               workc(i2+2,j1,k1) = aijk(2,i,j,kmkm)
            end do
            workc(imax+1,j1,k1) = 0.0d0
            workc(imax+2,j1,k1) = 0.0d0
         end do
         do j=jmax/2,jmax-1
            j1 = j+1
            jmjm = j-jmax
            do i=0,imax/2-1
               i2 = i+i
               workc(i2+1,j1,k1) = aijk(1,i,jmjm,kmkm)
               workc(i2+2,j1,k1) = aijk(2,i,jmjm,kmkm)
            end do
            workc(imax+1,j1,k1) = 0.0d0
            workc(imax+2,j1,k1) = 0.0d0
         end do
      end do
!$omp end parallel do
      
   end function aijk_to_workc
!-----------------------------------------------------------
end module flowfield
!===========================================================
module time_integration
   use,intrinsic :: iso_fortran_env
   use constants, only : nu,imax,jmax,kmax,ip,jp,kp
   implicit none
   private
   public timeintegration_initialize,rk4_for_nonlinear_ift_for_linear
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1) :: cl1ijk
   contains
   subroutine timeintegration_initialize(dt)
      use,intrinsic :: iso_fortran_env
      use flowfield, only : cn2ijk
      implicit none
      real(real64),intent(in) :: dt
      integer(int32) :: i,j,k


! constants for linear term calculation
      do k=-kmax/2,kmax/2-1
         do j=-jmax/2,jmax/2-1
            do i=0,imax/2-1
               cl1ijk(i,j,k) = exp(-0.50d0*dt*nu*cn2ijk(i,j,k))
            end do
         end do
      end do

! constants for nonlinear term calculation
!  #3
!  #4

      return
   end subroutine timeintegration_initialize
!-----------------------------------------------------------
   subroutine rk4_for_nonlinear_ift_for_linear(dt,f0,des_p_r2c,des_p_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      use flowfield,only : guarantee_solenoidal
      implicit none
      real(real64),intent(in) :: dt
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
      real(real64),parameter :: c16 = 1.0d0/6.0d0
      real(real64),parameter :: c13 = 1.0d0/3.0d0
      real(real64),dimension(imax*jmax*kmax*3),intent(inout) :: f0
      real(real64),dimension(imax*jmax*kmax*3) :: f1,f2,f3

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
      
      call guarantee_solenoidal(f0)

      return
   end subroutine rk4_for_nonlinear_ift_for_linear
!-----------------------------------------------------------
   subroutine intg_nonlinear(voijk,dvoijk,des_p_r2c,des_p_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      use fft_mkl, only : fft_execute_forward,fft_execute_backward
      use flowfield, only : aijk_to_workc,workc_to_aijk,p_length,vorticity_to_velocity
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(in)  :: voijk
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(out) :: dvoijk
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3) :: veijk
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
      real(real64),dimension(2,0:ip/2-1,-jp/2:jp/2-1,-kp/2:kp/2-1) :: ht_p
      real(real64),dimension(2,0:ip/2-1,-jp/2:jp/2-1,-kp/2:kp/2-1,3) :: vxvijk
      real(real64),dimension(ip+2,jp,kp) :: workc
      real(real64),dimension(ip,jp,kp,3) :: vorticity,velocity,vexvo
      integer(int32) :: i,j,k,n

! vorticity
      do n=1,3
!$omp parallel workshare
         ht_p(:,:,:,:) = 0.0d0
!$omp end parallel workshare
!$omp parallel workshare
         ht_p(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1)  &
&           =voijk(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,n)
!$omp end parallel workshare
         workc = aijk_to_workc(ip,jp,kp,ht_p)
         call fft_execute_backward(des_p_c2r,p_length,workc,vorticity(:,:,:,n))
      end do

! velocity
      call vorticity_to_velocity(voijk,veijk)
      do n=1,3
!$omp parallel workshare
         ht_p(:,:,:,:) = 0.0d0
!$omp end parallel workshare
!$omp parallel workshare
         ht_p(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1)  &
&           =veijk(1:2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,n)
!$omp end parallel workshare
         workc = aijk_to_workc(ip,jp,kp,ht_p)
         call fft_execute_backward(des_p_c2r,p_length,workc,velocity(:,:,:,n))
      end do

! calc. (velocity)X(vorticity)
      Block
      real(real64) :: u1,u2,u3,o1,o2,o3
!$omp parallel do shared(vorticity,velocity,vexvo),private(i,j,k,o1,o2,o3,u1,u2,u3)
      do k=1,kp
         do j=1,jp
            do i=1,ip
               o1 = vorticity(i,j,k,1)
               o2 = vorticity(i,j,k,2)
               o3 = vorticity(i,j,k,3)
               u1 = velocity (i,j,k,1)
               u2 = velocity (i,j,k,2)
               u3 = velocity (i,j,k,3)
               vexvo(i,j,k,1) = u2*o3-u3*o2
               vexvo(i,j,k,2) = u3*o1-u1*o3
               vexvo(i,j,k,3) = u1*o2-u2*o1
            end do
         end do
      end do
!$omp end parallel do
      end Block

! physical space to fourier
      do n=1,3
         call fft_execute_forward(des_p_r2c,p_length,vexvo(:,:,:,n),workc)
         vxvijk(:,:,:,:,n) = workc_to_aijk(ip,jp,kp,workc)
      end do

      Block
      real(real64) :: di,dj,dk
      real(real64),dimension(2,3) :: vxv
!$omp parallel do shared(vxvijk,dvoijk),private(i,j,k,di,dj,dk,vxv)
      do k=-kmax/2,kmax/2-1
         dk = dble(k)
         do j=-jmax/2,jmax/2-1
            dj = dble(j)
            do i=0,imax/2-1
               di = dble(i)
               vxv(1:2,1:3) = vxvijk(1:2,i,j,k,1:3)
               dvoijk(1,i,j,k,1) = -dj*vxv(2,3)+dk*vxv(2,2)
               dvoijk(2,i,j,k,1) =  dj*vxv(1,3)-dk*vxv(1,2)
               dvoijk(1,i,j,k,2) = -dk*vxv(2,1)+di*vxv(2,3)
               dvoijk(2,i,j,k,2) =  dk*vxv(1,1)-di*vxv(1,3)
               dvoijk(1,i,j,k,3) = -di*vxv(2,2)+dj*vxv(2,1)
               dvoijk(2,i,j,k,3) =  di*vxv(1,2)-dj*vxv(1,1)
            end do
         end do
      end do
!$omp end parallel do
      end Block

!!$omp parallel workshare
!      dvoijk(:,:,:,:,:) = -dvoijk(:,:,:,:,:)
!!$omp end parallel workshare
    
      return

   end subroutine intg_nonlinear
!-----------------------------------------------------------
   subroutine intg_linear(f)
      use,intrinsic :: iso_fortran_env
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1,-kmax/2:kmax/2-1,3),intent(inout) :: f
!      integer(int32) :: i,j,k

!$omp parallel workshare
      f(1,:,:,:,1) = f(1,:,:,:,1)*cl1ijk(:,:,:)
      f(2,:,:,:,1) = f(2,:,:,:,1)*cl1ijk(:,:,:)
      f(1,:,:,:,2) = f(1,:,:,:,2)*cl1ijk(:,:,:)
      f(2,:,:,:,2) = f(2,:,:,:,2)*cl1ijk(:,:,:)
      f(1,:,:,:,3) = f(1,:,:,:,3)*cl1ijk(:,:,:)
      f(2,:,:,:,3) = f(2,:,:,:,3)*cl1ijk(:,:,:)
!$omp end parallel workshare

      return

   end subroutine intg_linear
!-----------------------------------------------------------
end module time_integration
!===========================================================
program main
   use,intrinsic :: iso_fortran_env
!$ use omp_lib
   use constants
   use arrays
   use mkl_dfti
   use fft_mkl, only : fft_initialize,fft_finalize
   use flowfield, only : flowfield_initialize,wrtd
   use time_integration, only : timeintegration_initialize,rk4_for_nonlinear_ift_for_linear
   implicit none
   type(dfti_descriptor),pointer :: des_n_r2c, des_n_c2r ! normal    (length=N)
   type(dfti_descriptor),pointer :: des_p_r2c, des_p_c2r ! padding   (length=N*3/2)
   integer(int32) :: n
   real(real64) :: time
!$ integer(int32) :: nt
!$ real(real64) :: t1,t2,t3

!$omp parallel
!$ nt = omp_get_num_threads()
!$omp end parallel
!$ print *,'main) omp no. of threads :',nt

! pre-process
   print *,'main) initialization start'
   print *,'main) --fft--'
   print *,'main) for N to N transform'
   call fft_initialize(imax,jmax,kmax,des_n_r2c,des_n_c2r)
   print *,'main) for N*3/2 to N*3/2 transform'
   call fft_initialize(ip,jp,kp,des_p_r2c,des_p_c2r)
   print *,'main) --flow field--'
   open(20,file='./output/velocity_error.txt',form='formatted',status='replace')
   open(21,file='./output/vorticity_error.txt',form='formatted',status='replace')
   open(22,file='./output/tede.txt',form='formatted',status='replace')
   write(22,*) '# time - Energy - Dissipation rate - Enstrophy '
   call flowfield_initialize(voijk,des_n_r2c)
   print *,'main) --time integration--'
   call timeintegration_initialize(dt)
   print *,'main) initialization end'
!
   print *,'main) output flow field at t=0'
   call wrtd(0,voijk,des_n_c2r,des_p_r2c,des_p_c2r)

! main sequence
!$   t3 = 0.0d0
   do n=1,nmax
      time = dble(n-1)*dt
#ifdef debug
      print *,'loop',n
#endif
!$    t1 = omp_get_wtime()
      call rk4_for_nonlinear_ift_for_linear(dt,voijk,des_p_r2c,des_p_c2r)
!$    t2 = omp_get_wtime()
!$      t3 = t3+t2-t1
      if(mod(n,nout).eq.0) then
         call wrtd(n,voijk,des_n_c2r,des_p_r2c,des_p_c2r)
      end if
   end do

! post-process
!$   print *,'time',t3
   call fft_finalize(des_n_r2c,des_n_c2r)
   call fft_finalize(des_p_r2c,des_p_c2r)
   close(20)
   close(21)

   stop
end program main
