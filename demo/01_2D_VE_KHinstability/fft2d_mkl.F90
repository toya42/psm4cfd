module fft2d_mkl
   use,intrinsic :: iso_fortran_env
   use mkl_dfti
   implicit none
   private
   public fft2d_initialize,fft2d_finalize,fft2d_execute_forward,fft2d_execute_backward
   contains
!-----------------------------------------------------------
   subroutine fft2d_initialize(imax,jmax,des_r2c,des_c2r)
      use,intrinsic :: iso_fortran_env
      implicit none
      integer(int32),intent(in) :: imax,jmax
      type(dfti_descriptor),pointer :: des_r2c,des_c2r
#ifdef 32bit_integer
      integer(int32),dimension(12) :: dftistatus
      integer(int32),dimension(2) :: length
      integer(int32),dimension(3) :: stride_i,stride_o
      integer(int32) :: i
#elif 64bit_integer
      integer(int64),dimension(12) :: dftistatus
      integer(int64),dimension(2) :: length
      integer(int64),dimension(3) :: stride_i,stride_o
      integer(int64) :: i
#endif
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

      print *,'fft2d_initialize) imax   =',imax
      print *,'fft2d_initialize) jmax   =',jmax
!#ifdef debug
      do i=1,12
         if(dftistatus(i).ne.0) then
            print *,'fft2d_initialize) dft setting error:',i,dftistatus(i)
            if(.not. dftierrorclass(dftistatus(i),dfti_no_error)) then
               print *,'Error: ',dftierrormessage(dftistatus(i))
            end if
            stop
         end if
      end do
      print *,'fft2d_initialize) dft settings are completed'
!#endif

      return
   end subroutine fft2d_initialize
!-----------------------------------------------------------
   subroutine fft2d_finalize(des_r2c,des_c2r)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_r2c,des_c2r
#ifdef 32bit_integer
      integer(int32),dimension(2) :: dftistatus
      integer(int32) :: i
#elif 64bit_integer
      integer(int64),dimension(2) :: dftistatus
      integer(int64) :: i
#endif

      dftistatus(1)=DftiFreeDescriptor(Des_r2c)
      dftistatus(2)=DftiFreeDescriptor(Des_c2r)

!#ifdef debug
      do i=1,2
         if(dftistatus(i).ne.0) then
            print *,'fft2d_finalize) dft finalization error:',i,dftistatus(i)
            if(.not. dftierrorclass(dftistatus(i),dfti_no_error)) then
               print *,'Error: ',dftierrormessage(dftistatus(i))
            end if
            stop
         end if
      end do
      print *,'fft2d_finalize) dft finalization is completed'
!#endif
      return
   end subroutine fft2d_finalize
!-----------------------------------------------------------
   subroutine fft2d_execute_forward(des_r2c,length,real2d,complex2d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_r2c
      integer,dimension(2),intent(in) :: length
      real(real64),dimension(length(1)),intent(in)  :: real2d
      real(real64),dimension(length(2)),intent(out) :: complex2d
#ifdef 32bit_integer
      integer(int32) :: dftistatus
#elif 64bit_integer
      integer(int64) :: dftistatus
#endif

      dftistatus = dfticomputeforward(des_r2c,real2d,complex2d)
#ifdef debug
      if(dftistatus.ne.0) then
         print *,'fft2d_execute_forward) dft finalization error:',dftistatus
         if(.not. dftierrorclass(dftistatus,dfti_no_error)) then
            print *,'Error: ',dftierrormessage(dftistatus)
         end if
         stop
      end if
#endif
      return
   end subroutine fft2d_execute_forward
!-----------------------------------------------------------
   subroutine fft2d_execute_backward(des_c2r,length,complex2d,real2d)
      use,intrinsic :: iso_fortran_env
      implicit none
      type(dfti_descriptor),pointer :: des_c2r
      integer,dimension(2),intent(in) :: length
      real(real64),dimension(length(2)),intent(in)  :: complex2d
      real(real64),dimension(length(1)),intent(out) :: real2d
#ifdef 32bit_integer
      integer(int32) :: dftistatus
#elif 64bit_integer
      integer(int64) :: dftistatus
#endif

      dftistatus = dfticomputebackward(des_c2r,complex2d,real2d)
#ifdef debug
      if(dftistatus.ne.0) then
         print *,'fft2d_execute_forward) dft finalization error:',dftistatus
         if(.not. dftierrorclass(dftistatus,dfti_no_error)) then
            print *,'Error: ',dftierrormessage(dftistatus)
         end if
         stop
      end if
#endif
      return
   end subroutine fft2d_execute_backward
!-----------------------------------------------------------
end module fft2d_mkl

