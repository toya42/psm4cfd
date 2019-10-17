!===========================================================
! governing equation : 2d vorticity equation
! bondary condition : periodic
! numerical method
!     spatial discretization : spectral method
!        trial and test function : Fourier series
!     non-linear term : pseudo spectral method
!        de-aliasing : 3/2-rule (padding)
!     time integration : Runge-Kutta 4 (nonlinear term) 
!                        integrating factor (linear term)
!===========================================================
module constants
   use,intrinsic :: iso_fortran_env
   implicit none
   real(real64),parameter :: pi = 4.0d0*atan(1d0)
   real(real64),parameter :: nu = 2.0d-6

   integer(int32),parameter :: imax = 4*2**6
   integer(int32),parameter :: jmax = imax
   integer(int32),parameter :: ip = imax*3/2
   integer(int32),parameter :: jp = ip
   real(real64),parameter :: dt = 5.0d-3
   real(real64),parameter :: tmax = 20.0d0
!   integer(int32),parameter :: nmax = 4000   ! tmax/dt
!   integer(int32),parameter :: nout = 50
   integer(int32),parameter :: nmax = 10
   integer(int32),parameter :: nout = 10
end module constants
!===========================================================
module arrays
   use,intrinsic :: iso_fortran_env
   use constants
   implicit none
   real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: aij
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
   subroutine fft_initialize(imax,jmax,des_r2c,des_c2r)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      type(dfti_descriptor),pointer :: des_r2c,des_c2r
      integer(int32),dimension(12) :: dftistatus
      integer(int32),dimension(2) :: length
      integer(int32),dimension(3) :: stride_i,stride_o
      integer(int32) :: i
      real(real64) :: ni

      ni = 1.0d0/dble(imax*jmax)
      length(1) = imax
      length(2) = jmax
      stride_i(1) = 0
      stride_i(2) = 1
      stride_i(3) = imax
      stride_o(1) = 0
      stride_o(2) = 1
      stride_o(3) = imax/2+1
!r2c
      dftistatus( 1) = dfticreatedescriptor(des_r2c,dfti_double,dfti_real,2,length)
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
!#ifdef debug
      do i=1,12
         if(dftistatus(i).ne.0) then
            print *,'fft_initialize) dft setting error:',i,dftistatus(i)
            stop
         end if
      end do
      print *,'fft_initialize) dft settings are completed'
!#endif

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

!#ifdef debug
      do i=1,2
         if(dftistatus(i).ne.0) then
            print *,'fft_finalize) dft finalization error:',i
            stop
         end if
      end do
      print *,'fft_finalize) dft finalization is completed'
!#endif
      return
   end subroutine fft_finalize
!-----------------------------------------------------------
   subroutine fft_execute_forward(des_r2c,length,real2d,complex2d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_r2c
      integer(int32),dimension(2),intent(in) :: length
      real(real64),dimension(length(1)),intent(in)  :: real2d
      real(real64),dimension(length(2)),intent(out) :: complex2d
      integer(int32) :: dftistatus

      dftistatus = dfticomputeforward(des_r2c,real2d,complex2d)
#ifdef debug
      if(dftistatus.ne.0) then
         print *,'fft_execute_forward) dft finalization error:',dftistatus
         stop
      end if
#endif
      return
   end subroutine fft_execute_forward
!-----------------------------------------------------------
   subroutine fft_execute_backward(des_c2r,length,complex2d,real2d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_c2r
      integer(int32),dimension(2),intent(in) :: length
      real(real64),dimension(length(2)),intent(in)  :: complex2d
      real(real64),dimension(length(1)),intent(out) :: real2d
      integer(int32) :: dftistatus

      dftistatus = dfticomputebackward(des_c2r,complex2d,real2d)
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
   use constants, only : pi,imax,jmax
   use fft_mkl,only : fft_execute_forward,fft_execute_backward
   implicit none
   private
   public flowfield_initialize,wrtd
   public aij_to_workc,workc_to_aij
   public cn1ij
   real(real64),dimension(:),allocatable :: xi,yj
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1) :: cn1ij,cn2ij
   integer(int32),dimension(2) :: n_length
   contains
!-----------------------------------------------------------
   subroutine flowfield_initialize(aij,des_n_r2c)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: aij
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: bij
      type(dfti_descriptor),pointer :: des_n_r2c
      real(real64),dimension(imax,jmax) :: zeta
!      real(real64),dimension(imax,jmax) :: psi
      real(real64),dimension(imax+2,jmax) :: workc
      integer(int32) :: i,j
!      real(real64) :: x,x1,x2,y,y1,y2,sigma,sgm2i,c
      real(real64) :: x,x1,x2,y
      
! constants for aij<->bij and  nonlinear term calculation
!  #1 and #2
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               cn2ij(i,j) = -dble(i*i+j*j)
         end do
      end do
      cn2ij(0,0) = 1.0d0
      cn1ij = 1.0d0/cn2ij


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

      call fft_execute_forward(des_n_r2c,n_length,zeta,workc)
      bij = workc_to_aij(imax,jmax,workc)

!$omp parallel workshare
      aij(1,:,:) = bij(1,:,:)*cn1ij(:,:)
      aij(2,:,:) = bij(2,:,:)*cn1ij(:,:)
!$omp end parallel workshare

      return
   end subroutine flowfield_initialize
!-----------------------------------------------------------
   subroutine wrtd(loop,aij,des_n_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      integer(int32),intent(in) :: loop
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: bij
      type(dfti_descriptor),pointer :: des_n_c2r
      real(real64),dimension(imax,jmax) :: zeta

      integer(int32) :: i,j
      real(real64),dimension(imax+2,jmax) :: workc
      real(real64),dimension(:,:,:,:),allocatable :: buf
      character(len=50) :: filename
      real(real64) :: dammy

      if(loop.eq.0) then
         allocate(buf(imax,jmax,1,3))
         do j=1,jmax
            do i=1,imax
               buf(i,j,1,1) = xi(i)
               buf(i,j,1,2) = yj(j)
               buf(i,j,1,3) = 0.0d0
            end do
         end do
         open(10,file='output/grid.xyz',form='unformatted',access='stream',status='replace')
         i=imax
         j=jmax
         write(10) i,j,1
         write(10) buf
         close(10)
         deallocate(buf)
         deallocate(xi,yj)
      end if
!$omp parallel workshare
      bij(1,:,:) = aij(1,:,:)*cn2ij(:,:)
      bij(2,:,:) = aij(2,:,:)*cn2ij(:,:)
!$omp end parallel workshare
      workc = aij_to_workc(imax,jmax,bij)
      call fft_execute_backward(des_n_c2r,n_length,workc,zeta)

!! function file
!      allocate(buf(imax,jmax,1,1))
!      do j=1,jmax
!         do i=1,imax
!            buf(i,j,1,1)=zeta(i,j)
!         end do
!      end do
!      write(filename,'("output/zeta_",i5.5,".fun")') loop
!      open(10,file=filename,form='unformatted',access='stream',status='replace')
!      i=imax
!      j=jmax
!      write(10) i,j,1,1
!      write(10) buf
!      close(10)
!      deallocate(buf)

! q file
      allocate(buf(imax,jmax,1,5))
      do j=1,jmax
         do i=1,imax
            buf(i,j,1,1)=zeta(i,j)
            buf(i,j,1,2)=0.0d0
            buf(i,j,1,3)=0.0d0
            buf(i,j,1,4)=0.0d0
            buf(i,j,1,5)=zeta(i,j)
         end do
      end do
      write(filename,'("output/zeta_",i5.5,".q")') loop
      open(11,file=filename,form='unformatted',access='stream',status='replace')
      i=imax
      j=jmax
      write(11) i,j,1
      write(11) dammy,dammy,dammy,dammy
      write(11) buf
      close(11)
      deallocate(buf)

      print *,'wrtd) flowfield output,loop=',loop

      return
   end subroutine wrtd
!-----------------------------------------------------------
!   pure function workc_to_aij(imax,jmax,workc) result(aij)
   function workc_to_aij(imax,jmax,workc) result(aij)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),dimension(imax+2,jmax),intent(in) :: workc
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: aij
      integer(int32) :: i,j,i2,j1,jmjm

!$omp parallel do shared(aij,workc),private(i,j,j1,i2)
      do j=0,jmax/2-1
         j1 = j+1
         do i=0,imax/2-1
            i2 = i+i
            aij(1,i,j) = workc(i2+1,j1)
            aij(2,i,j) = workc(i2+2,j1)
         end do
      end do
!$omp end parallel do
!$omp parallel do shared(aij,workc),private(i,j,j1,jmjm,i2)
      do j=jmax/2,jmax-1
         j1 = j+1
         jmjm = j-jmax
         do i=0,imax/2-1
            i2 = i+i
            aij(1,i,jmjm) = workc(i2+1,j1) 
            aij(2,i,jmjm) = workc(i2+2,j1) 
         end do
      end do
!$omp end parallel do
      
   end function workc_to_aij
!-----------------------------------------------------------
!   pure function aij_to_workc(imax,jmax,aij) result(workc)
   function aij_to_workc(imax,jmax,aij) result(workc)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      real(real64),dimension(imax+2,jmax) :: workc
      integer(int32) :: i,j,i2,j1,jmjm

!$omp parallel do shared(workc,aij),private(i,j,j1,i2)
      do j=0,jmax/2-1
         j1 = j+1
         do i=0,imax/2-1
            i2 = i+i
            workc(i2+1,j1) = aij(1,i,j)
            workc(i2+2,j1) = aij(2,i,j)
         end do
         workc(imax+1,j1) = 0.0d0
         workc(imax+2,j1) = 0.0d0
      end do
!$omp end parallel do
!$omp parallel do shared(workc,aij),private(i,j,j1,i2,jmjm)
      do j=jmax/2,jmax-1
         j1 = j+1
         jmjm = j-jmax
         do i=0,imax/2-1
            i2 = i+i
            workc(i2+1,j1) = aij(1,i,jmjm)
            workc(i2+2,j1) = aij(2,i,jmjm)
         end do
         workc(imax+1,j1) = 0.0d0
         workc(imax+2,j1) = 0.0d0
      end do
!$omp end parallel do
      
   end function aij_to_workc
!-----------------------------------------------------------
end module flowfield
!===========================================================
module time_integration
   use,intrinsic :: iso_fortran_env
   use constants, only : nu,imax,jmax,ip,jp
   implicit none
   private
   public timeintegration_initialize,rk4_for_nonlinear_ift_for_linear
   real(real64),dimension(0:imax/2-1,-jmax/2:jmax/2-1) :: cl1ij,cn3ij,cn4ij
   integer(int32),dimension(2) :: p_length
   contains
   subroutine timeintegration_initialize(dt)
      use,intrinsic :: iso_fortran_env
      use flowfield, only : cn1ij
      implicit none
      real(real64),intent(in) :: dt
      integer(int32) :: i,j

      p_length(1) = ip*jp
      p_length(2) = (ip/2+1)*2*jp   ! (ip+2)*jp

! constants for linear term calculation
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            cl1ij(i,j) = dexp(-0.50d0*dt*nu*dble(i*i+j*j))
         end do
      end do

! constants for nonlinear term calculation
!  #3
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               cn3ij(i,j) = dble(i*i-j*j)*cn1ij(i,j)
         end do
      end do
!  #4
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
               cn4ij(i,j) = dble(i*j)*cn1ij(i,j)
         end do
      end do

      return
   end subroutine timeintegration_initialize
!-----------------------------------------------------------
   subroutine rk4_for_nonlinear_ift_for_linear(dt,f0,des_p_r2c,des_p_c2r)
      use,intrinsic :: iso_fortran_env
      use mkl_dfti
      implicit none
      real(real64),intent(in) :: dt
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
      real(real64),parameter :: c16 = 1.0d0/6.0d0
      real(real64),parameter :: c13 = 1.0d0/3.0d0
      real(real64),dimension(imax*jmax),intent(inout) :: f0
      real(real64),dimension(imax*jmax) :: f1,f2,f3

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
      use fft_mkl, only : fft_execute_forward,fft_execute_backward
      use flowfield, only : aij_to_workc,workc_to_aij
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: f
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: df
      type(dfti_descriptor),pointer :: des_p_r2c,des_p_c2r
      real(real64),dimension(2,0:ip/2-1,-jp/2:jp/2-1) :: ht_p
      real(real64),dimension(ip+2,jp) :: workc
      real(real64)   ,dimension(0:ip-1,0:jp-1)   :: u,v,uv,v2mu2
      integer(int32) :: i,j
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
      call fft_execute_backward(des_p_c2r,p_length,workc,u)

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
      call fft_execute_backward(des_p_c2r,p_length,workc,v)

! calc. Fij
!$omp parallel workshare
      uv(:,:)    = u(:,:)*v(:,:)
!$omp end parallel workshare
      call fft_execute_forward(des_p_r2c,p_length,uv,workc)
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
      call fft_execute_forward(des_p_r2c,p_length,v2mu2,workc)
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
   subroutine intg_linear(w)
      use,intrinsic :: iso_fortran_env
      implicit none
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(inout) :: w
      integer(int32) :: i,j

!$omp parallel do shared(w,cl1ij),private(i,j)
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            w(1,i,j) = w(1,i,j)*cl1ij(i,j)
            w(2,i,j) = w(2,i,j)*cl1ij(i,j)
         end do
      end do
!$omp end parallel do
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
   call fft_initialize(imax,jmax,des_n_r2c,des_n_c2r)
   print *,'main) for N*3/2 to N*3/2 transform'
   call fft_initialize(ip,jp,des_p_r2c,des_p_c2r)
!  initial flow field setting. zeta(vorticity) -> aij(spectral coefficients)
   print *,'main) --flow field--'
   call flowfield_initialize(aij,des_n_r2c)
   print *,'main) --time integration--'
   call timeintegration_initialize(dt)
   print *,'main) initialization end'
!
   print *,'main) output flow field at t=0'
   call wrtd(0,aij,des_n_c2r)

! main sequence
!$   t3 = 0.0d0
   do n=1,nmax
      time = dble(n-1)*dt
#ifdef debug
      print *,'loop',n
#endif
!$    t1 = omp_get_wtime()
      call rk4_for_nonlinear_ift_for_linear(dt,aij,des_p_r2c,des_p_c2r)
!$    t2 = omp_get_wtime()
!$      t3 = t3+t2-t1
      if(mod(n,nout).eq.0) then
         call wrtd(n,aij,des_n_c2r)
      end if
   end do

! post-process
!$   print *,'time',t3
   call fft_finalize(des_n_r2c,des_n_c2r)
   call fft_finalize(des_p_r2c,des_p_c2r)

   stop
end program main
