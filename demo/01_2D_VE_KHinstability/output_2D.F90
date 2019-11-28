module output
   use,intrinsic :: iso_fortran_env
   use constants, only : imax,jmax
   implicit none
   private
   public wrtd!,energy_spectrum_2D
   contains
!-----------------------------------------------------------
   subroutine wrtd(loop,aij,des_n_c2r)
      use fft2d_mkl,only : fft2d_execute_forward,fft2d_execute_backward
      use mkl_dfti
      use sort_spectral_coefficients, only : aij_to_workc
      use arrays, only : xi,yj,n_length,cn2ij
      implicit none
      integer(int32),intent(in) :: loop
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      type(dfti_descriptor),pointer :: des_n_c2r
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: bij
      real(real64),dimension(imax,jmax) :: zeta
#if integertype==0
      integer(int32) :: i,j,k
#elif integertype==1
      integer(int64) :: i,j,k
#endif
      integer(int32) :: io,jo,ko,lo
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
         io=imax
         jo=jmax
         ko=1
         write(10) io,jo,ko
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
      call fft2d_execute_backward(des_n_c2r,n_length,workc,zeta)

! function file
!      allocate(buf(imax,jmax,1,1))
!      do j=1,jmax
!         do i=1,imax
!            buf(i,j,1,1)=zeta(i,j)
!         end do
!      end do
!      write(filename,'("output/zeta_",i5.5,".fun")') loop
!      open(10,file=filename,form='unformatted',access='stream',status='replace')
!      io=imax
!      jo=jmax
!      ko=1
!      lo=1
!      write(10) io,jo,ko,lo
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
            buf(i,j,1,5)=0.0d0
         end do
      end do
      write(filename,'("output/zeta_",i5.5,".q")') loop
      open(11,file=filename,form='unformatted',access='stream',status='replace')
      io=imax
      jo=jmax
      ko=1
      write(11) io,jo,ko
      write(11) dammy,dammy,dammy,dammy
      write(11) buf
      close(11)
      deallocate(buf)

      print *,'wrtd) flowfield output,loop=',loop

      return
   end subroutine wrtd
!-----------------------------------------------------------
end module output

