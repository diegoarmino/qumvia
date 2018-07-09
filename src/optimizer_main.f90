subroutine optimize(qva_nml,nqmatoms,at_numbers,qmcoords,xopt)
  use qvamod_common, only: qva_nml_type
  use opt_data_mod, only : dXold_vec, dXnew_vec, Xold_vec, Xnew_vec, Hold, &
                            lambda, max_opt_steps
  implicit none
! -----------------------------------------------------------------------------
  type(qva_nml_type), intent(in)                  :: qva_nml  ! QUMVIA namelist parameters
  integer,            intent(in)                  :: nqmatoms
  integer,            intent(in)                  :: at_numbers(nqmatoms)
  double precision,   intent(in)                  :: qmcoords(3,nqmatoms)
  double precision,   intent(out)                 :: xopt(3,nqmatoms)
! -----------------------------------------------------------------------------

  include "qvmbia_param.f"
! -----------------------------------------------------------------------------

  max_opt_steps=qva_nml%max_opt_steps
  nat=nqmatoms

  allocate( Xnew(3,nat) )
  allocate( dXnew(3,nat) )

  Xnew=qmcoords
  dX=0d0
  Xold=0d0
  select case (qva_nml%opt_type)
    case ('SD')
    case ('CG')
      allocate( vold(3,nat)   )
    case ('QN')
      allocate( Xold_vec(3*nat)   )
      allocate( Xnew_vec(3*nat)   )
      allocate( dXnew_vec(3*nat)  )
      allocate( dXold_vec(3*nat)  )
      allocate( Hold(3*nat,3*nat) )
    case default
      allocate( Xold_vec(3*nat)   )
      allocate( Xnew_vec(3*nat)   )
      allocate( dXnew_vec(3*nat)  )
      allocate( dXold_vec(3*nat)  )
      allocate( Hold(3*nat,3*nat) )
  end select

   converged = .FALSE.
   do i=1,max_opt_steps
!     Single point energy and gradient evaluation.
      stepno = i
      call SCF_in(escf,Xnew,clcoords,nclatoms,dXnew)
      call dft_get_qm_forces(dXnew)
      call dft_get_mm_forces(dxyzcl,dXnew)

!     Check optimization convergence
      call opt_check_convergence(nat,at_numbers)
      if (convergence == .TRUE.) then
        call print_geom()
        return
      end if

      select case (qva_nml%opt_type)
      case ('SD')
        call steepest_descent(nqmatoms,xopt,at_numbers)
      case ('CG')
        call conjugated_gradient()
      case ('QN')
        call quasi_newton()
      case default
        call steepest_descent(nqmatoms,xopt,at_numbers)
      end select
   end do

   return
end subroutine optimize
