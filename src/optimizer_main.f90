subroutine optimize(nqmatoms,qmcoords,xopt)
   implicit none
!  -----------------------------------------------------------------------------
   integer,          intent(in)                  :: nqmatoms
   double precision, intent(in)                  :: qmcoords(3,nqmatoms)
   double precision, intent(out)                 :: xopt(3,nqmatoms)
!  -----------------------------------------------------------------------------
   double precision             :: lambda
   integer,parameter   :: nclatoms=0
   integer             :: ncoords              ! 3*nqmatoms
   integer             :: i,j,k,r,s,m
   double precision    :: dxyzcl(3,nclatoms)   ! SCF MM force
   double precision    :: dxyzqm(3,nqmatoms)   ! SCF QM force
   double precision    :: gradient(3*nqmatoms)

   include "qumvia_param.f90"
!  -----------------------------------------------------------------------------

   x0 = qmcoords/a0

   do i=1,max_opt_steps
!     Single point energy and gradient evaluation.
      call SCF_in(escf,qmcoords,clcoords,nclatoms,dipxyz)
      call dft_get_qm_forces(dxyzqm)
      call dft_get_mm_forces(dxyzcl,dxyzqm)

!     Change gradient into vector format.
      do i=1,nqmatoms
         do j=1,3
            gradient(3*(i-1)+j)=dxyzqm(j,i)
         end do
      end do

!     Check convergence
      call check_opt_convergence(converged)
!     Move to new geometry
      if (.NOT. converged) then
         call opt_mover()
      else
         write(77,'(A)') "Geometry Converged!"
         call write_geom(xopt)
         exit
      end if
   end do

   return
end subroutine optimize

subroutine opt_mover(nqmatoms,qmcoords,gradient,Xmoved)
   implicit none
   integer,          intent(in)                  :: nqmatoms
   double precision, intent(in)                  :: qmcoords(3,nqmatoms)
   double precision, intent(in)                  :: gradient(3,nqmatoms)
   double precision, intent(out)                 :: xmoved(3,nqmatoms)

   select case (optimizer)
   case ('SD')
      call steepest_descent()
   case ('CG')
      call conjugated_gradient()
   case default
      call steepest_descent()
   end select

end subroutine opt_mover

subroutine check_opt_convergence(nqmatoms,gradient,maxg_thresh,rms_thresh,converged)
   implicit none
   integer,          intent(in)                  :: nqmatoms
   double precision, intent(in)                  :: gradient(3*nqmatoms)
   logical,          intent(out)                 :: converged

   converged=F

   rms = dnrm2(3*nqmatoms,gradient,1)
   igmax = idamax(3*nqmatoms,gradient,1)

   if(gradient(igmax) < maxg_thresh .AND. rms < rms_thresh) then
      converged = T
   end if

end subroutine check_opt_convergence
