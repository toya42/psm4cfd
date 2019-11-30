module output
   use,intrinsic :: iso_fortran_env
   use constants, only : imax,jmax,pi
   implicit none
   private
   public wrtd,energy_spectra
   contains
!-----------------------------------------------------------
   subroutine energy_spectra(loop,aij)
      use arrays, only : cn2ij
      implicit none
      integer(int32),intent(in) :: loop
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: bij
#if integertype==0
      integer(int32) :: i,j,k
#elif integertype==1
      integer(int64) :: i,j,k
#endif
      character(len=50) :: filename
      integer(int32) :: kmax
      real(real64),dimension(:),allocatable :: ek
      real(real64) :: dk

      kmax=int(sqrt(dble((imax/2)**2+(jmax/2)**2)))
      allocate(ek(0:kmax))

!$omp parallel workshare
      bij(1,:,:) = aij(1,:,:)*cn2ij(:,:)
      bij(2,:,:) = aij(2,:,:)*cn2ij(:,:)
!$omp end parallel workshare

      ek = 0.0d0
      do j=-jmax/2,jmax/2-1
         do i=0,imax/2-1
            dk = sqrt(dble(i*i+j*j))
            k = int(dk)
            if(dk.eq.0.0d0) cycle
            ek(k) = ek(k)+pi/dk*(bij(1,i,j)**2+bij(2,i,j)**2)
         end do
      end do
      write(filename,'("output/energy_spectra_",i5.5,".txt")') loop
      open(12,file=filename,form='formatted',status='replace')
      write(12,*) '#',loop
      do k=0,kmax
         if(ek(k).ne.0) then
            write(12,'(i4,E18.8e3)') k,ek(k)
         end if
      end do
      close(12)

      return
   end subroutine energy_spectra
!-----------------------------------------------------------
   subroutine wrtd(oflag,loop,aij,des_n_c2r)
! oflag
! 0: no output    
! 1: grid file(.xyz)
! 2: grid file(.xyz) + one   variable(zeta)     function files(.fun)
! 3: grid file(.xyz) + one   variable(zeta)     function files(.fun) + dammy Q file(.q)
! 4: grid file(.xyz) + three variable(zeta,u,v) function files(.fun) 
! 5: grid file(.xyz) + three variable(zeta,u,v) function files(.fun) + dammy Q file(.q)
! 6: grid file(.xyz) + Q files(rho=zeta, rho*u=zeta*u, rho*v=zeta*v,w=e=0) (for paraview)
      use fft2d_mkl,only : fft2d_execute_forward,fft2d_execute_backward
      use mkl_dfti
      use sort_spectral_coefficients, only : aij_to_workc
      use arrays, only : xi,yj,n_length,cn2ij
      use velocity, only : psi_to_velocity
      implicit none
      integer(int32),intent(in) :: oflag,loop
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
      type(dfti_descriptor),pointer :: des_n_c2r
      real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1) :: bij
      real(real64),dimension(imax,jmax) :: zeta,vel1,vel2
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

      if(oflag==0) return

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
         
         if(oflag==3 .or. oflag==5) then
            allocate(buf(imax,jmax,1,5))
            buf(:,:,:,:)=0.0d0
            write(filename,'("output/dammy.q")')
            open(11,file=filename,form='unformatted',access='stream',status='replace')
            io=imax
            jo=jmax
            ko=1
            dammy=0.0d0
            write(11) io,jo,ko
            write(11) dammy,dammy,dammy,dammy
            write(11) buf
            close(11)
            deallocate(buf)
         end if
      end if

      if(oflag==1) return

      call psi_to_velocity(des_n_c2r,aij,vel1,vel2)
      
!$omp parallel workshare
      bij(1,:,:) = aij(1,:,:)*cn2ij(:,:)
      bij(2,:,:) = aij(2,:,:)*cn2ij(:,:)
!$omp end parallel workshare
      workc = aij_to_workc(imax,jmax,bij)
      call fft2d_execute_backward(des_n_c2r,n_length,workc,zeta)


! function file
      if(oflag==2 .or. oflag==3) then
         allocate(buf(imax,jmax,1,1))
         do j=1,jmax
            do i=1,imax
               buf(i,j,1,1)=zeta(i,j)
            end do
         end do
         write(filename,'("output/zeta_",i5.5,".fun")') loop
         open(10,file=filename,form='unformatted',access='stream',status='replace')
         io=imax
         jo=jmax
         ko=1
         lo=1
         write(10) io,jo,ko,lo
         write(10) buf
         close(10)
         deallocate(buf)
         print *,'wrtd) flowfield output,loop=',loop
         return
      end if

      if(oflag==4 .or. oflag==5) then
         allocate(buf(imax,jmax,1,3))
         do j=1,jmax
            do i=1,imax
               buf(i,j,1,1)=zeta(i,j)
               buf(i,j,1,2)=vel1(i,j)
               buf(i,j,1,3)=vel2(i,j)
            end do
         end do
         write(filename,'("output/zeta_u_v_",i5.5,".fun")') loop
         open(10,file=filename,form='unformatted',access='stream',status='replace')
         io=imax
         jo=jmax
         ko=1
         lo=3
         write(10) io,jo,ko,lo
         write(10) buf
         close(10)
         deallocate(buf)
         print *,'wrtd) flowfield output,loop=',loop
         return
      end if

! Q file
      if(oflag==6) then
         allocate(buf(imax,jmax,1,5))
         do j=1,jmax
            do i=1,imax
               buf(i,j,1,1)=zeta(i,j)
               buf(i,j,1,2)=zeta(i,j)*vel1(i,j)
               buf(i,j,1,3)=zeta(i,j)*vel2(i,j)
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
      end if
      
      if(oflag>6) then
         print *,'wrtd) oflag must be 0-6, but oflag=',oflag
         stop
      end if
   
   end subroutine wrtd
!-----------------------------------------------------------
end module output

