module flowfield
   use,intrinsic :: iso_fortran_env
   use constants, only : imax,jmax,pi
   implicit none
   private
   public flowfield_initialize_2D_VE_DP
! demo04

!-----------------------------------------------------------
   interface 
      module subroutine flowfield_initialize_2D_VE_DP(aij,des_n_r2c)
         use mkl_dfti
         implicit none
         real(real64),dimension(2,0:imax/2-1,-jmax/2:jmax/2-1),intent(out) :: aij
         type(dfti_descriptor),pointer :: des_n_r2c
      end subroutine
   end interface
!-----------------------------------------------------------
end module flowfield

