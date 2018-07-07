subroutine optimize(qva_nml,nqmatoms,at_numbers,qmcoords,xopt)
   use qvamod_common, only: qva_nml_type
   use opt_data_mod, only: max_opt_steps,nat,X0,Xold,dX
   implicit none
!  -----------------------------------------------------------------------------
   type(qva_nml_type), intent(in)                  :: qva_nml  ! QUMVIA namelist parameters
   integer,            intent(in)                  :: nqmatoms
   integer,            intent(in)                  :: at_numbers(nqmatoms)
   double precision,   intent(in)                  :: qmcoords(3,nqmatoms)
   double precision,   intent(out)                 :: xopt(3,nqmatoms)
!  -----------------------------------------------------------------------------

   include "qvmbia_param.f"
!  -----------------------------------------------------------------------------

   max_opt_steps=qva_nml%max_opt_steps
   nat=nqmatoms

   allocate( X0(3,nat) )
   allocate( dX(3,nat) )
   allocate( Xold(3,nat) )

   X0=qmcoords
   dX=0d0
   Xold=0d0

   do i=1,max_opt_steps
!     Single point energy and gradient evaluation.
      stepno = i
      call SCF_in(escf,X0,clcoords,nclatoms,dipxyz)
      call dft_get_qm_forces(dX)
      call dft_get_mm_forces(dxyzcl,dX)

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
