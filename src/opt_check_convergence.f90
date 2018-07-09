subroutine opt_check_convergence(nat,at_numbers)
   use opt_data_mod, only : dXnew, converged, qva_nml
   use qvamod_common, only: qva_nml_type
   implicit none

!  -----------------------------------------------------------------------------
   integer,            intent(in)     :: nat
   integer,            intent(in)     :: at_numbers(nat)
!  -----------------------------------------------------------------------------
   integer                            :: j,k
   integer                            :: igmax
   double precision                   :: rms
!  -----------------------------------------------------------------------------
   double precision,external          :: DNRM2
   double precision,external          :: IDAMAX
!  -----------------------------------------------------------------------------

! Change gradient into vector format.
  do j=1,nat
  do k=1,3
    dXnew_vec(3*(j-1)+k)   =  dXnew(k,j)
  end do
  end do

! Check convergence
  rms   = DNRM2(3*nat,dXnew_vec,1)
  rms   = rms/SQRT(3*DBLE(nat))
  igmax = idamax(3*nat,dXnew_vec,1)

  if(dX_vec(igmax) < qva_nml%maxg_thresh .AND. rms < qva_nml%rms_thresh) converged = .TRUE.

end subroutine conjugate_gradient
