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

  call opt_initialize(qva_nml)
  Xnew=qmcoords

  converged = .FALSE.
  do i=1,max_opt_steps

!   Single point energy and gradient evaluation.
    stepno = i
    call SCF_in(escf,Xnew,clcoords,nclatoms,dXnew)
    call dft_get_qm_forces(dXnew)
    call dft_get_mm_forces(dxyzcl,dXnew)

!   Check optimization convergence
    call opt_check_convergence(nat,at_numbers,converged)

!   If converged print optimized geometry and quit.
    if (converged == .TRUE.) then
      write(77,'(A)') "Geometry Converged!"
      do j=1,nat
         write(77,'(I3,3D18.8)') at_numbers(j),(Xnew(k,j),k=1,3)
      end do
      Xopt = Xnew
      return
    end if

!   Select optimization algorithm.
    select case (qva_nml%opt_type)
      case ('SD')
        call steepest_descent()
      case ('CG')
        call conjugated_gradient()
      case ('QN')
        call quasi_newton()
      case default
        call steepest_descent()
    end select

  end do

  return
end subroutine optimize
