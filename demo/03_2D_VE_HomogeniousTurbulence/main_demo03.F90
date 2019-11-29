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
program main
   use,intrinsic :: iso_fortran_env
!$ use omp_lib
   use constants
   use arrays, only : constant_parameters_initialize
   use mkl_dfti
   use fft2d_mkl, only : fft2d_initialize,fft2d_finalize
   use initial_flowfield, only : flowfield_initialize
   use output, only : wrtd
   use time_integration, only : rk4_for_nonlinear_ift_for_linear
   implicit none
   real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: aij
   type(dfti_descriptor),pointer :: des_n_r2c, des_n_c2r ! normal    (length=N)
   type(dfti_descriptor),pointer :: des_p_r2c, des_p_c2r ! padding   (length=N*3/2)
   integer(int32) :: n
   real(real64) :: time
!$ integer(int32) :: nt
!$ real(real64) :: t1,t2,t3

!$omp parallel
!$ nt = omp_get_num_threads()
!$omp end parallel
!$ print *,'main) OpenMP No. of threads :',nt

! pre-process
   print *,'main) pre-process start'
   print *,'main) --fft2d--'
   print *,'main) for N to N transform'
   call fft2d_initialize(imax,jmax,des_n_r2c,des_n_c2r)
   print *,'main) for N*3/2 to N*3/2 transform'
   call fft2d_initialize(ip,jp,des_p_r2c,des_p_c2r)
   print *,'main) --constant parameters--'
   call constant_parameters_initialize
   print *,'main) --flow field--'
   call flowfield_initialize(imax,jmax,aij,des_n_r2c)
   print *,'main) output flow field at t=0'
   n=0
   call wrtd(oflag,n,aij,des_n_c2r)
   print *,'main) pre-process end'

   print *,'main) main sequence start'
! main sequence
!$   t3 = 0.0d0
   do n=1,nmax
      time = dble(n-1)*dt
#ifdef debug
      print *,'loop',n
#endif
!$    t1 = omp_get_wtime()
      call rk4_for_nonlinear_ift_for_linear(dt,ijmax,aij,des_p_r2c,des_p_c2r)
!$    t2 = omp_get_wtime()
!$      t3 = t3+t2-t1
      if(mod(n,nout).eq.0) then
         call wrtd(oflag,n,aij,des_n_c2r)
      end if
   end do
   print *,'main) main sequence end'

! post-process
   print *,'main) post-process start'
!$   print *,'computational time',t3
   call fft2d_finalize(des_n_r2c,des_n_c2r)
   call fft2d_finalize(des_p_r2c,des_p_c2r)
   print *,'main) post-process end'
   print *,'main) calculation is completed successfully'

   stop
end program main
