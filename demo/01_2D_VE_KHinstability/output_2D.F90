module output_2D
   use,intrinsic :: iso_fortran_env
   use constants, only : imax,jmax,pi
   implicit none
   private
   public wrtd_2D_VE_DP,energy_spectrum_2D

!-----------------------------------------------------------
   interface 
      module subroutine wrtd_2D_VE_DP(loop,aij,des_n_c2r)
         use mkl_dfti
         use sort_spectral_coefficients, only : aij_to_workc
         implicit none
         integer(int32),intent(in) :: loop
         real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
         type(dfti_descriptor),pointer :: des_n_c2r
      end subroutine wrtd_2D_VE_DP
      module subroutine energy_spectrum_2D(loop,aij,des_n_c2r)
         use mkl_dfti
         use sort_spectral_coefficients, only : aij_to_workc
         implicit none
         integer(int32),intent(in) :: loop
         real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(in) :: aij
         type(dfti_descriptor),pointer :: des_n_c2r
      end subroutine energy_spectrum_2D

   end interface
!-----------------------------------------------------------
end module output_2D

