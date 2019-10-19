module sort_spectral_coefficients
   use,intrinsic :: iso_fortran_env
   implicit none
   private
   public aij_to_workc,workc_to_aij
   contains
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
end module sort_spectral_coefficients
