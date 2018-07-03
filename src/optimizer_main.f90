subroutine optimize(qva_nml,nqmatoms,qmcoords,xopt)
   use opt_data_mod, only: max_opt_steps,nat,X0
   implicit none
!  -----------------------------------------------------------------------------
   type(qva_nml_type), intent(in)                  :: qva_nml  ! QUMVIA namelist parameters
   integer,            intent(in)                  :: nqmatoms
   double precision,   intent(in)                  :: qmcoords(3,nqmatoms)
   double precision,   intent(out)                 :: xopt(3,nqmatoms)
!  -----------------------------------------------------------------------------

   include "qumvia_param.f90"
!  -----------------------------------------------------------------------------

   max_opt_steps=qva_nml%max_opt_steps
   nat=nqmatoms
   X0=qmcoords

   select case (qva_nml%opt_type)
   case ('SD')
      call steepest_descent(xopt)
!   case ('CG')
!      call conjugated_gradient()
   case default
      call steepest_descent()
   end select

   return
end subroutine optimize
