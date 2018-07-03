subroutine steepest_descent(Xopt)
   use opt_data_mod, only : nat, X0, dX, Xk, lambda, max_opt_steps
   implicit none

!  -----------------------------------------------------------------------------
   double precision,   intent(out)    :: Xopt(3,nat)
!  -----------------------------------------------------------------------------
   logical                            :: converged
   integer                            :: nclatoms=0
   double precision                   :: escf,eold
   double precision                   :: Xold(3,nat)
   double precision                   :: dipxyz(3)
   double precision                   :: dxyzcl(3,nat)
   double precision                   :: clcoords(4,nat)
   double precision                   :: gradvec(3*nat)
!  -----------------------------------------------------------------------------
   double precision,parameter         :: maxg_thresh=4d-2
   double precision,parameter         :: rms_thresh=2d-2
!  -----------------------------------------------------------------------------

   converged=F
   lambda=2.5D-2
   do i=1,max_opt_steps
!     Single point energy and gradient evaluation.
      call SCF_in(escf,X0,clcoords,nclatoms,dipxyz)
      call dft_get_qm_forces(dX)
      call dft_get_mm_forces(dxyzcl,dX)

!     Check if energy rises and update lambda
      if (i > 1 .AND. escf < eold) lambda = lambda*1.2d0
      if (i > 1 .AND. escf > eold) then
         lambda = lambda*0.5d0
         X0 = Xold
         CYCLE
      end if

!     Change gradient into vector format.
      do i=1,nat
         do j=1,3
            gradvec(3*(i-1)+j)=dX(j,i)
         end do
      end do

!     Check convergence

      rms   = dnrm2(3*nat,gradvec,1)
      rms   = rms/sqrt(3*nat)
      igmax = idamax(3*nat,gradvec,1)

      if(gradient(igmax) < maxg_thresh .AND. rms < rms_thresh) converged = T

!     Move to new geometry
      if (.NOT. converged) then
         gnorm =   dnrm2(3*nat,gradvec,1)
         svec  = - gradvec/gnorm
         Xold  =   X0
         X0    =   X0+lambda*svec
         eold  =   escf
      else
         Xopt  =  X0
         write(77,'(A)') "Geometry Converged!"
         call write_geom(xopt)
         exit
      end if
   end do

end subroutine steepest_descent
